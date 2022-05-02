import numpy as np
try:
    from mpi4py import MPI
    _mpi_available = True
except:
    _mpi_available = False

from enum import IntEnum
import h5py
from hazel.codes import hazel_code
from hazel.exceptions import NumericalErrorHazel, NumericalErrorSIR
from tqdm import tqdm
import logging
import os
import time
import signal
import sys
import warnings

class tags(IntEnum):
    READY = 0
    DONE = 1
    EXIT = 2
    START = 3
    DONOTHING = 4
    TOTALTIME = 5

class Iterator(object):
    def __init__(self, use_mpi=False):
        
        # Initializations and preliminaries        
        self.use_mpi = use_mpi
        
        if (self.use_mpi):
            if (not _mpi_available):
                raise Exception("You need MPI and mpi4py installed in your system to use this option.")
            self.comm = MPI.COMM_WORLD   # get MPI communicator object
            self.size = self.comm.Get_size()        # total number of processes
            self.rank = self.comm.Get_rank()        # rank of this process
            self.status = MPI.Status()   # get MPI status object                        
            warnings.filterwarnings("ignore", category=RuntimeWarning)  # Ignore all RuntimeWarnings
        else:
            self.rank = 0            

        self.logger = logging.getLogger("iterator")
        self.logger.setLevel(logging.DEBUG)
        self.logger.handlers = []

        ch = logging.StreamHandler()
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

    def get_rank(self, n_workers=0):        
        if (self.use_mpi):
            if (n_workers >= self.size):
                raise Exception("Number of requested workers {0} is >= number number of available cores ({1})".format(n_workers, size))
        return self.rank
    
    def use_model(self, model=None):
        
        # Then broadcast        
        if (self.use_mpi):
            if (self.rank == 0):
                self.model = model
                
                if (self.model.verbose):
                    self.logger.info('Broadcasting models to all workers')

                self.comm.Barrier()
                self.comm.bcast(self.model, root=0)                
                self.comm.Barrier()                
            else:
                model = None

                self.comm.Barrier()
                self.model = self.comm.bcast(model, root=0)
                self.comm.Barrier()
                # Initialize pyhazel
                hazel_code._init()
                self.model.init_sir()

            # In MPI mode, reduce MPI verbosity to a minimum
            self.model.verbose = 0
                            
        else:
            self.model = model

    def nonmpi_work(self, start=0, end=None):
        """
        Do the synthesis/inversion for all pixels in the models

        Parameters
        ----------
        model : model
            Model to be synthesized

        Returns
        -------
        None
        """

        if (end is None):
            end = self.model.n_pixels
            
        self.model.open_output()

        if (self.model.working_mode == 'inversion'):
            for k, v in self.model.spectrum.items():
                v.open_observation()
                v.mask_handle.open()
        
        # Open all model files
        for k, v in self.model.atmospheres.items():
            v.model_handler.open()
        
        # Loop over all pixels doing the synthesis and saving the results
        for i in tqdm(range(start, end)):

            mask = []

            for k, v in self.model.spectrum.items():
                v.read_mask(pixel=i)
                mask.append(v.mask)

            # Do the inversion only if the mask is set to 1
            if all(x == 1 for x in mask):

                for k, v in self.model.atmospheres.items():
                    if (v.type == 'photosphere'):
                        args, ff, vmac = v.model_handler.read(pixel=i)
                        v.set_parameters(args, ff, vmac)
                    else:
                        args, ff = v.model_handler.read(pixel=i)
                        v.set_parameters(args, ff)
                    
                    v.init_reference()

                if (self.model.n_randomization > 1):
                    randomize = True
                else:
                    randomize = False

                for loop in range(self.model.n_randomization):

                    if (self.model.working_mode == 'inversion'):
                        # Read current spectrum and stray light
                        for k, v in self.model.spectrum.items():
                            v.read_observation(pixel=i)
                        self.model.invert(randomize=randomize, randomization_ind=loop)                        
                    else:
                        self.model.synthesize()
                        self.model.flatten_parameters_to_reference(cycle=0)

                    self.model.output_handler.write(self.model, pixel=i, randomization=loop)
        
        if (self.model.working_mode == 'inversion'):
            for k, v in self.model.spectrum.items():
                v.close_observation()

        self.model.close_output()
                                            

    def mpi_parent_work(self, start=0, end=None):
        """
        MPI parent work

        Parameters
        ----------
        start : int
            Initial pixel, by default 0
        end : int
            Final pixel, by default the number of pixels in the observations

        Returns
        -------
        None
        """
        if (end is None):
            end = self.model.n_pixels

        self.model.open_output()

        if (self.model.working_mode == 'inversion'):
            for k, v in self.model.spectrum.items():
                v.open_observation()
                v.mask_handle.open()
        
        # Loop over all pixels doing the synthesis/inversion and saving the results
        task_index = start
        num_workers = self.size - 1
        closed_workers = 0
        self.last_received = 0
        self.last_sent = 0

        self.time = [None] * num_workers

        # Open all model files
        for k, v in self.model.atmospheres.items():
            v.model_handler.open()

        self.logger.info("Starting calculation with {0} workers".format(num_workers))

        self.elapsed = 0.0
        self.avg_elapsed = 0.0
        
        with tqdm(initial=task_index, total=end, ncols=140) as pbar:
            while (closed_workers < num_workers):
                data_received = self.comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=self.status)
                source = self.status.Get_source()
                tag = self.status.Get_tag()

                self.workers = num_workers - closed_workers
                
                if tag == tags.READY:
                    # Worker is ready, send a task
                    if (task_index < end):

                        data_to_send = {'index': task_index}

                        mask = []

                        for k, v in self.model.spectrum.items():
                            v.read_mask(pixel=task_index)
                            mask.append(v.mask)

                        # Do the inversion only if the mask is set to 1
                        if all(x == 1 for x in mask):
                        
                            if (self.model.working_mode == 'inversion'):

                                # Read current model atmosphere
                                for k, v in self.model.atmospheres.items():
                                    if (v.type == 'photosphere'):                                        
                                        args, ff, vmac = v.model_handler.read(pixel=task_index)
                                        v.set_parameters(args, ff, vmac)
                                    else:
                                        args, ff = v.model_handler.read(pixel=task_index)
                                        v.set_parameters(args, ff)
                                    v.init_reference()
                                    data_to_send[k] = [v.reference, v.parameters, v.stray_profile]

                                # Read current spectrum and stray light
                                for k, v in self.model.spectrum.items():
                                    v.read_observation(pixel=task_index)
                                    
                                    data_to_send[k] = [v.obs, v.noise, v.los, v.boundary, v.mu]
                            else:
                                for k, v in self.model.atmospheres.items():
                                    if (v.type == 'photosphere'):
                                        args, ff, vmac = v.model_handler.read(pixel=task_index)
                                        data_to_send[k] = [args, ff, vmac]
                                    else:
                                        args, ff = v.model_handler.read(pixel=task_index)
                                        data_to_send[k] = [args, ff]

                            self.comm.send(data_to_send, dest=source, tag=tags.START)

                        else:
                            self.comm.send(data_to_send, dest=source, tag=tags.DONOTHING)
                        
                        task_index += 1
                        pbar.update(1)
                        self.last_sent = '{0} to {1}'.format(task_index, source)

                        if (self.elapsed == 'Numerical error'):
                            pbar.set_postfix(sent=self.last_sent, received=self.last_received, workers=self.workers, elapsed='Numerical error')
                        else:
                            pbar.set_postfix(sent=self.last_sent, received=self.last_received, workers=self.workers, elapsed='{0:6.3f} s <{1:6.3f} s>'.format(self.elapsed, self.avg_elapsed))
                    else:
                        self.comm.send(None, dest=source, tag=tags.EXIT)
                elif tag == tags.DONE:
                    index = data_received['index']

                    if (data_received['error'] == 0):

                    # Loop over randomizations
                        for loop in range(self.model.n_randomization):                        
                            label = 'randomization_{0}'.format(loop)

                        # Put the results passed through MPI into self.model for saving
                            for k, v in self.model.spectrum.items():                
                                v.stokes_cycle = data_received[label][k][0]
                                v.chi2_cycle = data_received[label][k][1]
                                v.bic_cycle = data_received[label][k][2]
                                v.aic_cycle = data_received[label][k][3]
                                if (v.interpolate_to_lr):
                                    v.stokes_lr_cycle = data_received[label][k][4]
                                
                            
                            for k, v in self.model.atmospheres.items():
                                v.reference_cycle, v.error_cycle, v.nodes_location_cycle = data_received[label][k]

                            self.model.output_handler.write(self.model, pixel=index, randomization=loop)
                        
                        self.last_received = '{0} from {1}'.format(index, source)
                        self.elapsed = data_received['elapsed']
                        self.avg_elapsed = self.avg_elapsed + (self.elapsed - self.avg_elapsed) / (index + 1)
                        pbar.set_postfix(sent=self.last_sent, received=self.last_received, workers=self.workers, elapsed='{0:6.3f} s <{1:6.3f} s>'.format(self.elapsed, self.avg_elapsed))

                    else:

                        self.last_received = '{0} from {1}'.format(index, source)
                        self.elapsed = 'Numerical error'
                        pbar.set_postfix(sent=self.last_sent, received=self.last_received, workers=self.workers, elapsed='Numerical error')
                    
                    # self.logger.info('Received {0}->{1}'.format(index, source))

                elif tag == tags.DONOTHING:
                    pass

                elif tag == tags.TOTALTIME:
                    index = data_received['index']                                        
                    self.time[index - 1] = data_received['total_time']

                elif tag == tags.EXIT:                    
                    closed_workers += 1
                    pbar.set_postfix(sent=self.last_sent, received=self.last_received, workers='{0} finished'.format(source), elapsed='{0:6.3f} s <{1:6.3f} s>'.format(self.elapsed, self.avg_elapsed))

        self.logger.info('\nSummary')        
        total_time = sum(self.time)
        for i in range(num_workers):
            self.logger.info(f'Worker {i+1:03d} - Total t={self.time[i]:6.3f} - pct: {100*self.time[i]/total_time:6.3f}% - optimal: {1.0/num_workers*100:6.3f}%')

    def mpi_workers_work(self):
        """
        MPI workers work

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        self.time = 0
        
        while True:
            self.comm.send(None, dest=0, tag=tags.READY)
            data_received = self.comm.recv(source=0, tag=MPI.ANY_TAG, status=self.status)

            tag = self.status.Get_tag()
            
            if tag == tags.START:
                task_index = data_received['index']

                data_to_send = {'index': task_index}

                t0 = time.time()
                
                if (self.model.working_mode == 'inversion'):
                    
                    for k, v in self.model.spectrum.items():
                        v.obs, v.noise, v.los, v.boundary, v.mu = data_received[k]

                    for k, v in self.model.atmospheres.items():
                        v.reference, v.parameters, v.stray_profile = data_received[k]

                    
                    # Check for randomization
                    if (self.model.n_randomization > 1):
                        randomize = True
                    else:
                        randomize = False                                            

                    # Loop over randomizations
                    for loop in range(self.model.n_randomization):

                        label = 'randomization_{0}'.format(loop)
                        data_to_send[label] = {}
                        
                        # Try to do the inversion
                        try:
                            self.model.invert(randomize=randomize, randomization_ind=loop)
                            data_to_send['error'] = 0
                            for k, v in self.model.spectrum.items():
                                if (v.interpolate_to_lr):
                                    data_to_send[label][k] = [self.model.spectrum[k].stokes_cycle, self.model.spectrum[k].chi2_cycle, self.model.spectrum[k].bic_cycle, self.model.spectrum[k].aic_cycle, self.model.spectrum[k].stokes_lr_cycle]
                                else:
                                    data_to_send[label][k] = [self.model.spectrum[k].stokes_cycle, self.model.spectrum[k].chi2_cycle, self.model.spectrum[k].bic_cycle, self.model.spectrum[k].aic_cycle]

                            for k, v in self.model.atmospheres.items():
                                data_to_send[label][k] = [v.reference_cycle, v.error_cycle, v.nodes_location_cycle]

                        # If a numerical problem appeared, send the error code to the parent
                        except NumericalErrorHazel:                            
                            data_to_send['error'] = 1
                            for k, v in self.model.spectrum.items():
                                data_to_send[label][k] = None

                            for k, v in self.model.atmospheres.items():
                                data_to_send[label][k] = None
                        
                        except NumericalErrorSIR:
                            data_to_send['error'] = 2
                            for k, v in self.model.spectrum.items():
                                data_to_send[label][k] = None

                            for k, v in self.model.atmospheres.items():
                                data_to_send[label][k] = None

                        except:
                            data_to_send['error'] = 3
                            for k, v in self.model.spectrum.items():
                                data_to_send[label][k] = None

                            for k, v in self.model.atmospheres.items():
                                data_to_send[label][k] = None
                        
                else:
                    for k, v in self.model.atmospheres.items():                    
                        if (v.type == 'photosphere'):
                            v.set_parameters(data_received[k][0], data_received[k][1], data_received[k][2])
                        else:
                            v.set_parameters(data_received[k][0], data_received[k][1])
                        
                    label = 'randomization_0'
                    data_to_send[label] = {}

                    try:
                        self.model.synthesize()
                        self.model.flatten_parameters_to_reference(cycle=0)
                        data_to_send['error'] = 0
                    except NumericalErrorHazel:                            
                        data_to_send['error'] = 1
                    except NumericalErrorSIR:                            
                        data_to_send['error'] = 2
                    except:
                        data_to_send['error'] = 3
                    
                    
                    for k, v in self.model.spectrum.items():
                        if (v.interpolate_to_lr):
                            data_to_send[label][k] = [self.model.spectrum[k].stokes_cycle, None, None, None, self.model.spectrum[k].stokes_lr_cycle]
                        else:
                            data_to_send[label][k] = [self.model.spectrum[k].stokes_cycle, None, None, None]

                    for k, v in self.model.atmospheres.items():
                        data_to_send[label][k] = [v.reference_cycle, v.error_cycle, None]
                                        
                t1 = time.time()
                data_to_send['elapsed'] = t1 - t0

                # Measure total computing time per node
                self.time += t1 - t0

                self.comm.send(data_to_send, dest=0, tag=tags.DONE)
            elif tag == tags.DONOTHING:
                self.comm.send({}, dest=0, tag=tags.DONOTHING)
            elif tag == tags.EXIT:
                data_to_send = {'index': self.rank, 'total_time': self.time}
                self.comm.send(data_to_send, dest=0, tag=tags.TOTALTIME)
                break

        self.comm.send(None, dest=0, tag=tags.EXIT)           

    def cleanup(self, *args):
        if (self.rank == 0):
            self.logger.info('Closing output file before exiting with exception')
            self.model.close_output()
            self.logger.info('Closed')
            sys.exit(0)

    def run_all_pixels(self, start=0, end=None):
        """
        Run synthesis/inversion for all pixels
        
        Parameters
        ----------
        start : int, optional
            Initial pixel, by default 0

        end : int
            Final pixel, by default the number of pixels in the observations
        """

        # Capture SIGTERM and SIGINT from the main process. When the process being run with mpiexec 
        # is killed, mpiexec 
        # During a normal python execution, hitting ^C sends a SIGINT signal to python, which 
        # translates it into a KeyboardInterrupt exception (https://docs.python.org/3.7/library/signal.html).
        # But when hitting ^C during an mpiexec execution, it is the mpiexec process which receives the 
        # SIGINT signal, and instead of propagating it to its children processes, it sends to its 
        # children processes a SIGTERM signal (https://www.open-mpi.org/doc/current/man1/mpirun.1.php).
        # It thus seems that python doesn't react similarly to SIGINT and SIGTERM signals.
        for sig in (signal.SIGINT, signal.SIGTERM):
            signal.signal(sig, self.cleanup)

        try:
            if (self.use_mpi):
                if (self.rank == 0):
                    self.mpi_parent_work(start, end)
                    self.logger.info('Finished.')
                else:
                    self.mpi_workers_work()
            else:
                self.nonmpi_work(start, end)

            if (self.rank == 0):
                try:
                    os.remove('lte.grid')
                except OSError:
                    pass
        except (KeyboardInterrupt, SystemExit):
            pass
