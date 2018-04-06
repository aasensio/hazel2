import numpy as np
from mpi4py import MPI
from enum import IntEnum
import h5py
from hazel.codes import hazel_code
from tqdm import tqdm
# from ipdb import set_trace as stop

class tags(IntEnum):
    READY = 0
    DONE = 1
    EXIT = 2
    START = 3

class iterator(object):
    def __init__(self, use_mpi=False):
        
        # Initializations and preliminaries
        self.use_mpi = use_mpi
        
        if (self.use_mpi):
            self.comm = MPI.COMM_WORLD   # get MPI communicator object
            self.size = self.comm.size        # total number of processes
            self.rank = self.comm.rank        # rank of this process
            self.status = MPI.Status()   # get MPI status object                        
        else:
            self.rank = 0            

    def get_rank(self, n_slaves=0):        
        if (self.use_mpi):
            if (n_slaves >= self.size):
                raise Exception("Number of requested slaves {0} is >= number number of available cores ({1})".format(n_slaves, size))        
        return self.rank
    
    def use_model(self, model=None):
        
        # Then broadcast        
        if (self.use_mpi):
            if (self.rank == 0):
                self.model = model

                if (self.model.verbose):
                    print ('Broadcasting models to all slaves')

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
                            
        else:
            self.model = model

    def nonmpi_work(self):
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

        self.model.open_output()
        
        # Open all model files
        for k, v in self.model.atmospheres.items():
            v.model_handler.open()
        
        # Loop over all pixels doing the synthesis and saving the results
        for i in tqdm(range(self.model.n_pixels)):
            for k, v in self.model.atmospheres.items():                
                args = v.model_handler.read(pixel=i)                
                v.set_parameters(args)

            if (self.model.working_mode == 'inversion'):
                # Read current spectrum and stray light
                for k, v in self.model.spectrum.items():
                    v.read_observation(pixel=i)
                    v.read_straylight(pixel=i)
                self.model.invert()
            else:
                self.model.synthesize()

            self.model.output_handler.write(self.model, pixel=i)

        self.model.close_output()
                                            

    def mpi_master_work(self):
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

        self.model.open_output()

        if (self.model.working_mode == 'inversion'):
            for k, v in self.model.spectrum.items():
                v.open_observation()
        
        # Loop over all pixels doing the synthesis/inversion and saving the results
        task_index = 0
        num_workers = self.size - 1
        closed_workers = 0
        self.last_received = 0
        self.last_sent = 0

        # Open all model files
        for k, v in self.model.atmospheres.items():
            v.model_handler.open()

        print("Starting calculation with {0} workers".format(num_workers), flush=True)
        
        with tqdm(total=self.model.n_pixels, ncols=140) as pbar:
            while (closed_workers < num_workers):
                data_received = self.comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=self.status)
                source = self.status.Get_source()
                tag = self.status.Get_tag()
                
                if tag == tags.READY:
                    # Worker is ready, send a task
                    if (task_index < self.model.n_pixels):

                        data_to_send = {'index': task_index}

                        if (self.model.working_mode == 'inversion'):
                            # Read current spectrum and stray light
                            for k, v in self.model.spectrum.items():
                                v.read_observation(pixel=task_index)
                                v.read_straylight(pixel=task_index)
                                data_to_send[k] = [v.obs, v.noise, v.stray, v.los, v.boundary]
                        else:
                            for k, v in self.model.atmospheres.items():
                                m = v.model_handler.read(pixel=task_index)
                                data_to_send[k] = m

                        self.comm.send(data_to_send, dest=source, tag=tags.START)

                        task_index += 1
                        pbar.update(1)
                        self.last_sent = '{0}->{1}'.format(task_index, source)
                        pbar.set_postfix(sent=self.last_sent, received=self.last_received)
                    else:
                        self.comm.send(None, dest=source, tag=tags.EXIT)
                elif tag == tags.DONE:
                    index = data_received['index']

                    # Put the results passed through MPI into self.model for saving
                    for k, v in self.model.spectrum.items():                
                        v.stokes = data_received[k]

                    self.model.output_handler.write(self.model, pixel=index)                    
                                    
                    # for k, v in self.model.spectrum.items():                    
                        # db[k][index,:,:] = data_received[k]                    
                    
                    self.last_received = '{0}->{1}'.format(index, source)
                    pbar.set_postfix(sent=self.last_sent, received=self.last_received)

                elif tag == tags.EXIT:                    
                    closed_workers += 1

        self.model.close_output()

    def mpi_slave_work(self):
        
        while True:
            self.comm.send(None, dest=0, tag=tags.READY)
            data_received = self.comm.recv(source=0, tag=MPI.ANY_TAG, status=self.status)

            tag = self.status.Get_tag()
            
            if tag == tags.START:                                
                task_index = data_received['index']

                
                
                if (self.model.working_mode == 'inversion'):
                    
                    for k, v in self.model.spectrum.items():
                        v.obs, v.noise, v.stray, v.los, v.boundary = data_received[k]

                    self.model.invert()
                else:
                    for k, v in self.model.atmospheres.items():                    
                        v.set_parameters(data_received[k])
                        
                    self.model.synthesize()

                data_to_send = {'index': task_index}

                for k, v in self.model.spectrum.items():
                    data_to_send[k] = self.model.spectrum[k].stokes
                    
                self.comm.send(data_to_send, dest=0, tag=tags.DONE)
            elif tag == tags.EXIT:
                break

        self.comm.send(None, dest=0, tag=tags.EXIT)           

    def run_all_pixels(self):
        if (self.use_mpi):
            if (self.rank == 0):
                self.mpi_master_work()
            else:
                self.mpi_slave_work()
        else:
            self.nonmpi_work()