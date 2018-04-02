import sys, os, random
import numpy as np
import i0Allen
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import pyhazel
from ipdb import set_trace as stop
import pickle
import os.path

import matplotlib.pyplot as pl
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class FloatSlider(QSlider):
    def __init__(self, parent, decimals=3, *args, **kargs):
        super(FloatSlider, self).__init__(parent, *args, **kargs)
        self._multi = 10 ** decimals
        self.setMinimum(self.minimum())
        self.setMaximum(self.maximum())

    def value(self):
        return float(super(FloatSlider, self).value()) / self._multi

    def setMinimum(self, value):
        return super(FloatSlider, self).setMinimum(value * self._multi)

    def setMaximum(self, value):
        return super(FloatSlider, self).setMaximum(value * self._multi)

    def setValue(self, value):
        super(FloatSlider, self).setValue(int(value * self._multi))

class ExtendedQLabel(QLabel):
 
    def __init(self, parent):
        QLabel.__init__(self, parent)
 
    def mouseReleaseEvent(self, ev):
        self.emit(SIGNAL('clicked()'))


class AppForm(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.setWindowTitle('hazel')

        self.resize(1200,200)

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()
        
        # self.on_draw()

    def saveConfig(self):
        d = {'synModeInput' : self.synModeInput, 'nSlabsInput' : self.nSlabsInput, 'B1Input' : self.B1Input,
            'B2Input' : self.B2Input, 'hInput' : self.hInput, 'tau1Input' : self.tau1Input,
            'tau2Input' : self.tau2Input, 'boundaryInput'  : self.boundaryInput,
            'transInput' : self.transInput, 'atomicPolInput' : self.atomicPolInput, 'anglesInput' : self.anglesInput,
            'lambdaAxisInput' : self.lambdaAxisInput, 'nLambdaInput' : self.nLambdaInput,
            'dopplerWidthInput' : self.dopplerWidthInput, 'dopplerWidth2Input' : self.dopplerWidth2Input, 'dampingInput' : self.dampingInput,
            'dopplerVelocityInput' : self.dopplerVelocityInput, 'dopplerVelocity2Input' : self.dopplerVelocity2Input, 'ffInput' : self.ffInput,
            'betaInput' : self.betaInput, 'beta2Input' : self.beta2Input, 'nbarInput' : self.nbarInput, 'omegaInput' : self.omegaInput, 'obsFile' : self.obsFile,
            'normalization': 0}
        pickle.dump( d, open( "state.pickle", "wb" ) )

    def loadConfig(self):
        if (os.path.exists("state.pickle")):
            try:
                d = pickle.load( open( "state.pickle", "rb" ) )
                self.synModeInput = d['synModeInput']
                self.nSlabsInput = d['nSlabsInput']
                self.B1Input = d['B1Input']
                self.B2Input = d['B2Input']
                self.hInput = d['hInput']
                self.tau1Input = d['tau1Input']
                self.tau2Input = d['tau2Input']
                self.boundaryInput  = d['boundaryInput']
                self.transInput = d['transInput']
                self.atomicPolInput = d['atomicPolInput']
                self.anglesInput = d['anglesInput']
                self.lambdaAxisInput = d['lambdaAxisInput']
                self.nLambdaInput = d['nLambdaInput']
                self.dopplerWidthInput = d['dopplerWidthInput']
                self.dopplerWidth2Input = d['dopplerWidth2Input']
                self.dampingInput = d['dampingInput']
                self.dopplerVelocityInput = d['dopplerVelocityInput']
                self.dopplerVelocity2Input = d['dopplerVelocity2Input']
                self.ffInput = d['ffInput']
                self.betaInput = d['betaInput']
                self.beta2Input = d['beta2Input']
                self.nbarInput = d['nbarInput']
                self.omegaInput = d['omegaInput']
                self.obsFile = d['obsFile']
                self.normalization = d['normalization']
            except:
                print('state.pickle from other python version')
                self.synModeInput = 5
                self.nSlabsInput = 1
                self.B1Input = np.asarray([3.0,80.0,41.0])
                self.B2Input = np.asarray([0.0,0.0,0.0])
                self.hInput = 3.e0
                self.tau1Input = 1.e0
                self.tau2Input = 0.e0
                self.boundaryInput  = np.asarray([0.0,0.0,0.0,0.0])
                self.transInput = 1
                self.atomicPolInput = 1
                self.anglesInput = np.asarray([90.0,0.0,90.0])
                self.lambdaAxisInput = np.asarray([-1.5e0,2.5e0])
                self.nLambdaInput = 150
                self.dopplerWidthInput = 6.e0
                self.dopplerWidth2Input = 0.e0
                self.dampingInput = 0.e0
                self.dopplerVelocityInput = 0.e0
                self.dopplerVelocity2Input = 0.e0
                self.ffInput = 0.e0
                self.betaInput = 1.0
                self.beta2Input = 1.0
                self.nbarInput = np.asarray([0.0,0.0,0.0,0.0])
                self.omegaInput = np.asarray([0.0,0.0,0.0,0.0])
                self.obsFile = ''
                self.normalization = 0
        else:
            self.synModeInput = 5
            self.nSlabsInput = 1
            self.B1Input = np.asarray([3.0,80.0,41.0])
            self.B2Input = np.asarray([0.0,0.0,0.0])
            self.hInput = 3.e0
            self.tau1Input = 1.e0
            self.tau2Input = 0.e0
            self.boundaryInput  = np.asarray([0.0,0.0,0.0,0.0])
            self.transInput = 1
            self.atomicPolInput = 1
            self.anglesInput = np.asarray([90.0,0.0,90.0])
            self.lambdaAxisInput = np.asarray([-1.5e0,2.5e0])
            self.nLambdaInput = 150
            self.dopplerWidthInput = 6.e0
            self.dopplerWidth2Input = 0.e0
            self.dampingInput = 0.e0
            self.dopplerVelocityInput = 0.e0
            self.dopplerVelocity2Input = 0.e0
            self.ffInput = 0.e0
            self.betaInput = 1.0
            self.beta2Input = 1.0
            self.nbarInput = np.asarray([0.0,0.0,0.0,0.0])
            self.omegaInput = np.asarray([0.0,0.0,0.0,0.0])
            self.obsFile = ''
            self.normalization = 0
        
    
#######################################################################
# EVENTS
#######################################################################
    def on_about(self):
        msg = """ Hazel:
        
         * A. Asensio Ramos
         * Instituto de Astrofisica de Canarias, Spain
        """
        QMessageBox.about(self, "About hazel", msg.strip())
    
    # def on_pick(self, event):
    #     # The event received here is of the type
    #     # matplotlib.backend_bases.PickEvent
    #     #
    #     # It carries lots of information, of which we're using
    #     # only a small amount here.
    #     # 
    #     box_points = event.artist.get_bbox().get_points()
    #     msg = "You've clicked on a bar with coords:\n %s" % box_points
        
    #     QMessageBox.information(self, "Click!", msg)
    
    # def on_draw(self):
    #     """ Redraws the figure
    #     """
    #     self.data = np.arange(10)
        
    #     x = np.arange(10)

    #     # clear the axes and redraw the plot anew
    #     #
    #     self.axes.clear()        
    #     # self.axes.grid(self.grid_cb.isChecked())
        
    #     self.axes.bar(
    #         left=x, 
    #         height=self.data, 
    #         width=self.slider.value() / 100.0, 
    #         align='center', 
    #         alpha=0.44,
    #         picker=5)
        
    #     self.canvas.draw()

    def redrawProfiles(self):
        self.status_text = "Computing"
        lambdaAxisInputA = np.linspace(self.lambdaAxisInput[0],self.lambdaAxisInput[1],self.nLambdaInput)
        [l, stokes, etaOutput, epsOutput] = pyhazel.synth(self.synModeInput, self.nSlabsInput, self.B1Input, self.B2Input, self.hInput, 
                        self.tau1Input, self.tau2Input, self.boundaryInput, self.transInput, self.atomicPolInput, self.anglesInput, 
                        self.nLambdaInput, lambdaAxisInputA, self.dopplerWidthInput, self.dopplerWidth2Input, self.dampingInput, 
                        self.dopplerVelocityInput, self.dopplerVelocity2Input, self.ffInput, self.betaInput, self.beta2Input, self.nbarInput, self.omegaInput,
                        self.normalization)
        
        for i in range(4):         
            self.axes[i].clear()
            self.axes[i].plot(l - 10829.0911,stokes[i,:])
            if (self.obsFile != ''):
                self.axes[i].plot(l - 10829.0911, self.obsStokes[:,i+1], 'r.')
            for item in ([self.axes[i].title, self.axes[i].xaxis.label, self.axes[i].yaxis.label] +
                self.axes[i].get_xticklabels() + self.axes[i].get_yticklabels()):
                item.set_fontsize(10)

        self.fig.canvas.draw()
        self.status_text = "OK"

    def loadObservation(self):
        self.obsFile = QFileDialog.getOpenFileName(self, 'Open file', '')
        if (self.obsFile != ''):
            self.obsStokes = np.loadtxt(str(self.obsFile))
            self.loadedFile.setText('Loaded: {0}'.format(self.obsFile))

    def resetObservation(self):
        self.obsFile = ''
        self.obsStokes = None
        self.loadedFile.setText('Loaded: none')

# B1
    def onSliderB1(self):
        self.sliderValueB1.setText('{0:6.1f}'.format(self.sliderB1.value()))
        self.B1Input[0] = self.sliderB1.value()
        self.redrawProfiles()    

    def onSliderLeftB1(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter lower limit:')
        if (ok):
            self.sliderB1.setMinimum(float(text))
            self.sliderLeftB1.setText(text)

    def onSliderRightB1(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter upper limit:')
        if (ok):
            self.sliderB1.setMaximum(float(text))
            self.sliderRightB1.setText(text)

    def onSliderValueB1(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter value B1:')
        if (ok):
            self.sliderB1.setValue(float(text))
            self.sliderValueB1.setText(text)
            self.B1Input[0] = self.sliderB1.value()
            self.redrawProfiles()

# thetaB1
    def onSliderthB1(self):
        self.sliderValuethB1.setText(str(self.sliderthB1.value()))
        self.B1Input[1] = self.sliderthB1.value()
        self.redrawProfiles()

    def onSliderValuethB1(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter value thB1:')
        if (ok):
            self.sliderthB1.setValue(float(text))
            self.sliderValuethB1.setText(text)
            self.B1Input[1] = self.sliderthB1.value()
            self.redrawProfiles()

# phiB1
    def onSliderphiB1(self):
        self.sliderValuephiB1.setText(str(self.sliderphiB1.value()))
        self.B1Input[2] = self.sliderphiB1.value()
        self.redrawProfiles()

    def onSliderValuephiB1(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter value phiB1:')
        if (ok):
            self.sliderphiB1.setValue(float(text))
            self.sliderValuephiB1.setText(text)

            self.B1Input[2] = self.sliderphiB1.value()
            self.redrawProfiles()

# vth1
    def onSlidervthB1(self):
        self.sliderValuevthB1.setText(str(self.slidervthB1.value()))
        self.dopplerWidthInput = self.slidervthB1.value()
        self.redrawProfiles()

    def onSliderValuevthB1(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter value vth1:')
        if (ok):
            self.slidervthB1.setValue(float(text))
            self.sliderValuevthB1.setText(text)
            self.dopplerWidthInput = self.slidervthB1.value()
            self.redrawProfiles()

# v1
    def onSlidervB1(self):
        self.sliderValuevB1.setText(str(self.slidervB1.value()))
        self.dopplerVelocityInput = self.slidervB1.value()
        self.redrawProfiles()

    def onSliderValuevB1(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter value vth1:')
        if (ok):
            self.slidervB1.setValue(float(text))
            self.sliderValuevB1.setText(text)
            self.dopplerVelocityInput = self.slidervB1.value()
            self.redrawProfiles()

# tau1
    def onSlidertau1(self):
        self.sliderValuetau1.setText(str(self.slidertau1.value()))
        self.tau1Input = self.slidertau1.value()
        self.redrawProfiles()

    def onSliderValuetau1(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter value tau1:')
        if (ok):
            self.slidertau1.setValue(float(text))
            self.sliderValuetau1.setText(text)
            self.tau1Input = self.slidertau1.value()
            self.redrawProfiles()

# B2
    def onSliderB2(self):
        self.sliderValueB2.setText(str(self.sliderB2.value()))
        self.B2Input[0] = self.sliderB2.value()
        self.redrawProfiles()

    def onSliderLeftB2(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter lower limit:')
        if (ok):
            self.sliderB2.setMinimum(float(text))
            self.sliderLeftB2.setText(text)

    def onSliderRightB2(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter upper limit:')
        if (ok):
            self.sliderB2.setMaximum(float(text))
            self.sliderRightB2.setText(text)

    def onSliderValueB2(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter value:')
        if (ok):
            self.sliderB2.setValue(float(text))
            self.sliderValueB2.setText(text)
            self.B2Input[0] = self.sliderB2.value()
            self.redrawProfiles()

# thetaB2
    def onSliderthB2(self):
        self.sliderValuethB2.setText(str(self.sliderthB2.value()))
        self.B2Input[1] = self.sliderB2.value()
        self.redrawProfiles()

    def onSliderValuethB2(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter value:')
        if (ok):
            self.sliderthB2.setValue(float(text))
            self.sliderValuethB2.setText(text)
            self.B2Input[1] = self.sliderthB2.value()
            self.redrawProfiles()

# phiB2
    def onSliderphiB2(self):
        self.sliderValuephiB2.setText(str(self.sliderphiB2.value()))
        self.B2Input[2] = self.sliderphiB2.value()
        self.redrawProfiles()

    def onSliderValuephiB2(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter value:')
        if (ok):
            self.sliderphiB2.setValue(float(text))
            self.sliderValuephiB2.setText(text)
            self.B2Input[2] = self.sliderphiB2.value()
            self.redrawProfiles()

# vth2
    def onSlidervthB2(self):
        self.sliderValuevthB2.setText(str(self.slidervthB2.value()))
        self.dopplerWidthInput2 = self.slidervthB2.value()
        self.redrawProfiles()

    def onSliderValuevthB2(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter value vth2:')
        if (ok):
            self.slidervthB2.setValue(float(text))
            self.sliderValuevthB2.setText(text)
            self.dopplerWidthInput2 = self.slidervthB2.value()
            self.redrawProfiles()

# v2
    def onSlidervB2(self):
        self.sliderValuevB2.setText(str(self.slidervB2.value()))
        self.dopplerVelocityInput2 = self.slidervB2.value()
        self.redrawProfiles()

    def onSliderValuevB2(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter value vth2:')
        if (ok):
            self.slidervB2.setValue(float(text))
            self.sliderValuevB2.setText(text)
            self.dopplerVelocityInput2 = self.slidervB2.value()
            self.redrawProfiles()

# tau2
    def onSlidertau2(self):
        self.sliderValuetau2.setText(str(self.slidertau2.value()))
        self.tau2Input = self.slidertau2.value()
        self.redrawProfiles()

    def onSliderValuetau2(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter value tau2:')
        if (ok):
            self.slidertau2.setValue(float(text))
            self.sliderValuetau2.setText(text)
            self.tau2Input = self.slidertau2.value()
            self.redrawProfiles()

# theta
    def onSliderTheta(self):
        self.sliderValuetheta.setText(str(self.slidertheta.value()))
        if (self.checkAllen.isChecked()):
            theta = float(self.sliderValuetheta.text())
            mu = np.cos(theta * np.pi / 180.0)
            i0 = i0Allen.i0Allen(10830.0, mu)
            self.I0.setText(str(i0))
        self.anglesInput[0] = self.slidertheta.value()
        self.redrawProfiles()

    def onSliderValueTheta(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter theta angle:')
        if (ok):
            self.slidertheta.setValue(float(text))
            self.sliderValuetheta.setText(str(text))
            self.anglesInput[0] = self.slidertheta.value()
            self.redrawProfiles()

# phi
    def onSliderPhi(self):
        self.sliderValuephi.setText(str(self.sliderphi.value()))
        self.anglesInput[1] = self.sliderphi.value()
        self.redrawProfiles()

    def onSliderValuePhi(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter phi angle:')
        if (ok):
            self.sliderphi.setValue(float(text))
            self.sliderValuephi.setText(str(text))
            self.anglesInput[1] = self.sliderphi.value()
            self.redrawProfiles()

# gamma
    def onSliderGamma(self):
        self.sliderValuegamma.setText(str(self.slidergamma.value()))
        self.anglesInput[2] = self.slidergamma.value()
        self.redrawProfiles()

    def onSliderValueGamma(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter gamma angle:')
        if (ok):
            self.slidergamma.setValue(float(text))
            self.sliderValuegamma.setText(str(text))
            self.anglesInput[2] = self.slidergamma.value()
            self.redrawProfiles()

# height
    def onSliderHeight(self):
        self.sliderValueheight.setText(str(self.sliderheight.value()))
        self.hInput = self.sliderheight.value()
        self.redrawProfiles()

    def onSliderValueHeight(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter height:')
        if (ok):
            self.sliderheight.setValue(float(text))
            self.sliderValueheight.setText(str(text))
            self.hInput = self.sliderheight.value()
            self.redrawProfiles()

# damping
    def onSliderDamp(self):
        self.sliderValueDamp.setText(str(self.sliderDamp.value()))
        self.dampingInput = self.sliderDamp.value()
        self.redrawProfiles()

    def onSliderValueDamp(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter damping:')
        if (ok):
            self.sliderDamp.setValue(float(text))
            self.sliderValueDamp.setText(str(text))
            self.dampingInput = self.sliderDamp.value()
            self.redrawProfiles()

# ff
    def onSliderff(self):
        self.sliderValuerff.setText(str(self.sliderff.value()))
        self.ffInput = self.slideff.value()
        self.redrawProfiles()

    def onSliderValueff(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter ff:')
        if (ok):
            self.sliderff.setValue(float(text))
            self.sliderValueff.setText(str(text))
            self.ffInput = self.slideff.value()
            self.redrawProfiles()

# S2/S1 slab 1
    def onSliderS2S1(self):        
        self.sliderValueS2S1.setText(str(self.sliderS2S1.value()))
        self.betaInput = self.sliderS2S1.value()
        self.redrawProfiles()

    def onSliderValueS2S1(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter S2/S1 slab 1:')
        if (ok):
            self.sliderS2S1.setValue(float(text))
            self.sliderValueS2S1.setText(str(text))
            self.betaInput = self.sliderS2S1.value()
            self.redrawProfiles()

# S2/S1 slab 2
    def onSliderS2S1B(self):        
        self.sliderValueS2S1B.setText(str(self.sliderS2S1B.value()))
        self.beta2Input = self.sliderS2S1B.value()
        self.redrawProfiles()

    def onSliderValueS2S1B(self):
        text, ok = QInputDialog.getText(self, 'Input Dialog', 'Enter S2/S1 slab 2:')
        if (ok):
            self.sliderS2S1B.setValue(float(text))
            self.sliderValueS2S1B.setText(str(text))
            self.beta2Input = self.sliderS2S1B.value()
            self.redrawProfiles()

# multiplet
    def onSlidermultiplet(self):     
        self.sliderValuermultiplet.setText(self.multiplets[self.slidermultiplet.value()-1])

# thin
    def onRadioThin(self):     
        self.synModeInput = 0
        self.redrawProfiles()

# exact
    def onRadioExact(self):     
        self.synModeInput = 5
        self.atomicPolInput = 1
        self.redrawProfiles()
# atompol
    def onRadioAtompol(self):
        self.synModeInput = 5    
        self.atomicPolInput = 0
        self.redrawProfiles()

# maxnormalization
    def onRadioMaxNorm(self):
        self.normalization = 0
        self.redrawProfiles()

# maxnormalization
    def onRadioPeakNorm(self):
        self.normalization = 1
        self.redrawProfiles()


    def onCheckAllen(self, state):
        if (state == 2):
            self.sliderValuertheta.setText(str(self.slidertheta.value()))
            if (self.checkAllen.isChecked()):
                theta = float(self.sliderValuertheta.text())
                mu = np.cos(theta * np.pi / 180.0)
                i0 = i0Allen.i0Allen(10830.0, mu)
                self.I0.setText(str(i0))
            self.anglesInput[0] = self.slidertheta.value()
            self.I0.setEnabled(False)
            self.Q0.setEnabled(False)
            self.U0.setEnabled(False)
            self.V0.setEnabled(False)
        else:
            self.I0.setEnabled(True)
            self.Q0.setEnabled(True)
            self.U0.setEnabled(True)
            self.V0.setEnabled(True)

# I0, Q0, U0, V0
    def onChangeI0(self):        
        self.boundaryInput[0] = float(self.I0.text())

    def onChangeQ0(self):     
        self.boundaryInput[1] = float(self.I0.text())

    def onChangeU0(self):     
        self.boundaryInput[2] = float(self.I0.text())

    def onChangeV0(self):     
        self.boundaryInput[3] = float(self.I0.text())

    def onChangeleftWave(self):
        self.lambdaAxisInput[0] = float(self.leftWave.text())
        self.redrawProfiles()

    def onChangerightWave(self):
        self.lambdaAxisInput[1] = float(self.rightWave.text())
        self.redrawProfiles()

    def onChangestepWave(self):
        self.nLambdaInput = int(self.stepWave.text())
        self.redrawProfiles()

    def onChangenSlabs(self, index):
        if (index == 0):
            self.nSlabsInput = 1
        if (index == 1):
            self.nSlabsInput = 2
        if (index == 2):
            self.nSlabsInput = -2
        self.redrawProfiles()

    def onChangeSlabs(self):
        pass
    
#######################################################################
# INITIALIZATION
#######################################################################
    def create_main_frame(self):
        self.main_frame = QWidget()

        self.fontSize = 9
        self.font = QFont("SansSerif", self.fontSize)

        QApplication.setFont(self.font)
        
# Hazel configuration
        self.loadConfig()
        pyhazel.init()

        self.multiplets = ['10830', '3888', '7065',' 5876']
        
        # Create the mpl Figure and FigCanvas objects. 
        # 5x4 inches, 100 dots-per-inch
        #

        self.minSizeSlider = 100
        self.maxSizeSlider = 200

        self.minLabelSize = 50
        self.maxLabelSize = 150

        self.dpi = 80
        # self.fig, self.ax = pl.subplots(ncols=2, nrows=2, figsize=(6,6))
        self.fig = Figure((8.0, 7.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        layout = QVBoxLayout()
        # self.canvas.setParent(self.main_frame)
        layout.addWidget(self.canvas)

        self.axes = [None]*4

        for i in range(4):
            self.axes[i] = self.fig.add_subplot(2,2,i+1)

        self.fig.tight_layout()

        
        # Bind the 'pick' event for clicking on one of the bars
        #
        # self.canvas.mpl_connect('pick_event', self.on_pick)
        
        # Create the navigation toolbar, tied to the canvas
        #
        # self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        
        # First component
        # 
        boxB1 = QGridLayout()        

        sliderLabel = QLabel('B [G]')
        self.sliderB1 = FloatSlider(Qt.Horizontal)
        self.sliderLeftB1 = ExtendedQLabel('0')
        self.sliderRightB1 = ExtendedQLabel('1000')
        self.sliderB1.setMinimum(0.0)
        self.sliderB1.setMaximum(1000.0)
        self.sliderB1.setValue(self.B1Input[0])
        self.sliderB1.setTracking(True)
        self.sliderB1.setTickPosition(QSlider.TicksBothSides)
        self.sliderB1.setMinimumWidth(self.minSizeSlider)
        self.sliderB1.setMaximumWidth(self.maxSizeSlider)
        self.sliderValueB1 = ExtendedQLabel('{0:6.1f}'.format(self.B1Input[0]))
        
        for l, w in enumerate([sliderLabel, self.sliderLeftB1, self.sliderB1, self.sliderRightB1, self.sliderValueB1]):
            boxB1.addWidget(w, 0, l)
            boxB1.setAlignment(w, Qt.AlignVCenter)
        self.sliderValueB1.setMinimumWidth(self.minLabelSize)
        self.sliderValueB1.setMaximumWidth(self.maxLabelSize)
        self.connect(self.sliderB1, SIGNAL('sliderReleased()'), self.onSliderB1)
        self.connect(self.sliderLeftB1, SIGNAL('clicked()'), self.onSliderLeftB1)
        self.connect(self.sliderRightB1, SIGNAL('clicked()'), self.onSliderRightB1)
        self.connect(self.sliderValueB1, SIGNAL('clicked()'), self.onSliderValueB1)
                

        sliderLabel = QLabel('thB [deg]')
        self.sliderthB1 = FloatSlider(Qt.Horizontal)
        self.sliderLeftthB1 = ExtendedQLabel('0')
        self.sliderRightthB1 = ExtendedQLabel('180')
        self.sliderthB1.setMinimum(0.0)
        self.sliderthB1.setMaximum(180.0)
        self.sliderthB1.setValue(self.B1Input[1])
        self.sliderthB1.setTracking(True)
        self.sliderthB1.setTickPosition(QSlider.TicksBothSides)
        self.sliderthB1.setMinimumWidth(self.minSizeSlider)
        self.sliderthB1.setMaximumWidth(self.maxSizeSlider)
        self.sliderValuethB1 = ExtendedQLabel('{0:6.1f}'.format(self.B1Input[1]))
        for l, w in enumerate([sliderLabel, self.sliderLeftthB1, self.sliderthB1, self.sliderRightthB1, self.sliderValuethB1]):
            boxB1.addWidget(w, 1, l)
            boxB1.setAlignment(w, Qt.AlignVCenter)
        self.sliderValuethB1.setMinimumWidth(self.minLabelSize)
        self.sliderValuethB1.setMaximumWidth(self.maxLabelSize)
        self.connect(self.sliderthB1, SIGNAL('sliderReleased()'), self.onSliderthB1)
        self.connect(self.sliderValuethB1, SIGNAL('clicked()'), self.onSliderValuethB1)

        sliderLabel = QLabel('phiB [deg]')
        self.sliderphiB1 = FloatSlider(Qt.Horizontal)
        self.sliderLeftphiB1 = ExtendedQLabel('-180')
        self.sliderRightphiB1 = ExtendedQLabel('180')
        self.sliderphiB1.setMinimum(-180)
        self.sliderphiB1.setMaximum(180)
        self.sliderphiB1.setValue(20)
        self.sliderphiB1.setTracking(True)
        self.sliderphiB1.setTickPosition(QSlider.TicksBothSides)
        self.sliderphiB1.setMinimumWidth(self.minSizeSlider)
        self.sliderphiB1.setMaximumWidth(self.maxSizeSlider)
        self.sliderValuephiB1 = ExtendedQLabel('{0:6.1f}'.format(self.B1Input[2]))
        for l, w in enumerate([sliderLabel, self.sliderLeftphiB1, self.sliderphiB1, self.sliderRightphiB1, self.sliderValuephiB1]):
            boxB1.addWidget(w, 2, l)
            boxB1.setAlignment(w, Qt.AlignVCenter)
        self.sliderValuephiB1.setMinimumWidth(self.minLabelSize)
        self.sliderValuephiB1.setMaximumWidth(self.maxLabelSize)
        self.connect(self.sliderphiB1, SIGNAL('sliderReleased()'), self.onSliderphiB1)
        self.connect(self.sliderValuephiB1, SIGNAL('clicked()'), self.onSliderValuephiB1)

        sliderLabel = QLabel('width [km/s]')
        self.slidervthB1 = FloatSlider(Qt.Horizontal)
        self.sliderLeftvthB1 = ExtendedQLabel('0.10')
        self.sliderRightvthB1 = ExtendedQLabel('20')
        self.slidervthB1.setMinimum(0.1)
        self.slidervthB1.setMaximum(20.0)
        self.slidervthB1.setValue(self.dopplerWidthInput)
        self.slidervthB1.setTracking(True)
        self.slidervthB1.setTickPosition(QSlider.TicksBothSides)
        self.slidervthB1.setMinimumWidth(self.minSizeSlider)
        self.slidervthB1.setMaximumWidth(self.maxSizeSlider)
        self.sliderValuevthB1 = ExtendedQLabel('{0:6.1f}'.format(self.dopplerWidthInput))
        for l, w in enumerate([sliderLabel, self.sliderLeftvthB1, self.slidervthB1, self.sliderRightvthB1, self.sliderValuevthB1]):
            boxB1.addWidget(w, 3, l)
            boxB1.setAlignment(w, Qt.AlignVCenter)
        self.sliderValuevthB1.setMinimumWidth(self.minLabelSize)
        self.sliderValuevthB1.setMaximumWidth(self.maxLabelSize)
        self.connect(self.slidervthB1, SIGNAL('sliderReleased()'), self.onSlidervthB1)
        self.connect(self.sliderValuevthB1, SIGNAL('clicked()'), self.onSliderValuevthB1)

        sliderLabel = QLabel('v [km/s]')
        self.slidervB1 = FloatSlider(Qt.Horizontal)
        self.sliderLeftvB1 = ExtendedQLabel('-10')
        self.sliderRightvB1 = ExtendedQLabel('10')
        self.slidervB1.setMinimum(-10)
        self.slidervB1.setMaximum(10)
        self.slidervB1.setValue(self.dopplerVelocityInput)
        self.slidervB1.setTracking(True)
        self.slidervB1.setTickPosition(QSlider.TicksBothSides)
        self.slidervB1.setMinimumWidth(self.minSizeSlider)
        self.slidervB1.setMaximumWidth(self.maxSizeSlider)
        self.sliderValuevB1 = ExtendedQLabel('{0:6.1f}'.format(self.dopplerVelocityInput))
        for l, w in enumerate([sliderLabel, self.sliderLeftvB1, self.slidervB1, self.sliderRightvB1, self.sliderValuevB1]):
            boxB1.addWidget(w, 4, l)
            boxB1.setAlignment(w, Qt.AlignVCenter)
        self.sliderValuevB1.setMinimumWidth(self.minLabelSize)
        self.sliderValuevB1.setMaximumWidth(self.maxLabelSize)
        self.connect(self.slidervB1, SIGNAL('sliderReleased()'), self.onSlidervB1)
        self.connect(self.sliderValuevB1, SIGNAL('clicked()'), self.onSliderValuevB1)

        sliderLabel = QLabel('tau')
        self.slidertau1 = FloatSlider(Qt.Horizontal)
        self.sliderLefttau1 = ExtendedQLabel('0.0')
        self.sliderRighttau1 = ExtendedQLabel('10')
        self.slidertau1.setMinimum(0.0)
        self.slidertau1.setMaximum(10.0)
        self.slidertau1.setValue(self.tau1Input)
        self.slidertau1.setTracking(True)
        self.slidertau1.setTickPosition(QSlider.TicksBothSides)
        self.slidertau1.setMinimumWidth(self.minSizeSlider)
        self.slidertau1.setMaximumWidth(self.maxSizeSlider)
        self.sliderValuetau1 = ExtendedQLabel('{0:6.1f}'.format(self.tau1Input))
        for l, w in enumerate([sliderLabel, self.sliderLefttau1, self.slidertau1, self.sliderRighttau1, self.sliderValuetau1]):
            boxB1.addWidget(w, 5, l)
            boxB1.setAlignment(w, Qt.AlignVCenter)
        self.sliderValuetau1.setMinimumWidth(self.minLabelSize)
        self.sliderValuetau1.setMaximumWidth(self.maxLabelSize)
        self.connect(self.slidertau1, SIGNAL('sliderReleased()'), self.onSlidertau1)
        self.connect(self.sliderValuetau1, SIGNAL('clicked()'), self.onSliderValuetau1)

        # Second component
        # 
        boxB2 = QGridLayout()
        sliderLabel = QLabel('B [G]')
        self.sliderB2 = FloatSlider(Qt.Horizontal)
        self.sliderLeftB2 = ExtendedQLabel('0')
        self.sliderRightB2 = ExtendedQLabel('1000')
        self.sliderB2.setMinimum(0.0)
        self.sliderB2.setMaximum(1000.0)
        self.sliderB2.setValue(self.B2Input[0])
        self.sliderB2.setTracking(True)
        self.sliderB2.setTickPosition(QSlider.TicksBothSides)
        self.sliderValueB2 = ExtendedQLabel('{0:6.1f}'.format(self.B2Input[0]))
        for l, w in enumerate([sliderLabel, self.sliderLeftB2, self.sliderB2 , self.sliderRightB2, self.sliderValueB2]):
            boxB2.addWidget(w, 0, l)
            boxB2.setAlignment(w, Qt.AlignVCenter)
        self.sliderValueB2.setMinimumWidth(self.minLabelSize)
        self.sliderValueB2.setMaximumWidth(self.maxLabelSize)
        self.connect(self.sliderB2, SIGNAL('sliderReleased()'), self.onSliderB2)
        self.connect(self.sliderLeftB2, SIGNAL('clicked()'), self.onSliderLeftB2)
        self.connect(self.sliderRightB2, SIGNAL('clicked()'), self.onSliderRightB2)
        self.connect(self.sliderValueB2, SIGNAL('clicked()'), self.onSliderValueB2)

        sliderLabel = QLabel('thB [deg]')
        self.sliderthB2 = FloatSlider(Qt.Horizontal)
        self.sliderLeftthB2 = ExtendedQLabel('0')
        self.sliderRightthB2 = ExtendedQLabel('180')
        self.sliderthB2.setMinimum(0)
        self.sliderthB2.setMaximum(180)
        self.sliderthB2.setValue(self.B2Input[1])
        self.sliderthB2.setTracking(True)
        self.sliderthB2.setTickPosition(QSlider.TicksBothSides)
        self.sliderthB2.setMinimumWidth(self.minSizeSlider)
        self.sliderthB2.setMaximumWidth(self.maxSizeSlider)
        self.sliderValuethB2 = ExtendedQLabel('{0:6.1f}'.format(self.B2Input[1]))
        for l, w in enumerate([sliderLabel, self.sliderLeftthB2, self.sliderthB2, self.sliderRightthB2, self.sliderValuethB2]):
            boxB2.addWidget(w, 1, l)
            boxB2.setAlignment(w, Qt.AlignVCenter)
        self.sliderValuethB2.setMinimumWidth(self.minLabelSize)
        self.sliderValuethB2.setMaximumWidth(self.maxLabelSize)
        self.connect(self.sliderthB2, SIGNAL('sliderReleased()'), self.onSliderthB2)
        self.connect(self.sliderValuethB2, SIGNAL('clicked()'), self.onSliderValuethB2)        

        sliderLabel = QLabel('phiB [deg]')
        self.sliderphiB2 = FloatSlider(Qt.Horizontal)
        self.sliderLeftPhiB2 = ExtendedQLabel('0')
        self.sliderRightPhiB2 = ExtendedQLabel('180')
        self.sliderphiB2.setMinimum(-180)
        self.sliderphiB2.setMaximum(180)
        self.sliderphiB2.setValue(self.B2Input[2])
        self.sliderphiB2.setTracking(True)
        self.sliderphiB2.setTickPosition(QSlider.TicksBothSides)
        self.sliderphiB2.setMinimumWidth(self.minSizeSlider)
        self.sliderphiB2.setMaximumWidth(self.maxSizeSlider)
        self.sliderValuephiB2 = ExtendedQLabel('{0:6.1f}'.format(self.B2Input[2]))
        for l, w in enumerate([sliderLabel, self.sliderLeftPhiB2, self.sliderphiB2, self.sliderRightPhiB2, self.sliderValuephiB2]):
            boxB2.addWidget(w, 2, l)
            boxB2.setAlignment(w, Qt.AlignVCenter)
        self.sliderValuephiB2.setMinimumWidth(self.minLabelSize)
        self.sliderValuephiB2.setMaximumWidth(self.maxLabelSize)
        self.connect(self.sliderphiB2, SIGNAL('sliderReleased()'), self.onSliderphiB2)
        self.connect(self.sliderValuephiB2, SIGNAL('clicked()'), self.onSliderValuephiB2)

        sliderLabel = QLabel('width [km/s]')
        self.slidervthB2 = FloatSlider(Qt.Horizontal)
        self.sliderLeftvthB2 = ExtendedQLabel('0')
        self.sliderRightvthB2 = ExtendedQLabel('180')
        self.slidervthB2.setMinimum(0.1)
        self.slidervthB2.setMaximum(20)
        self.slidervthB2.setValue(self.dopplerWidth2Input)
        self.slidervthB2.setTracking(True)
        self.slidervthB2.setTickPosition(QSlider.TicksBothSides)
        self.slidervthB2.setMinimumWidth(self.minSizeSlider)
        self.slidervthB2.setMaximumWidth(self.maxSizeSlider)
        self.sliderValuevthB2 = ExtendedQLabel('{0:6.1f}'.format(self.dopplerWidth2Input))
        for l, w in enumerate([sliderLabel, self.sliderLeftvthB2, self.slidervthB2, self.sliderRightvthB2, self.sliderValuevthB2]):
            boxB2.addWidget(w, 3, l)
            boxB2.setAlignment(w, Qt.AlignVCenter)
        self.sliderValuevthB2.setMinimumWidth(self.minLabelSize)
        self.sliderValuevthB2.setMaximumWidth(self.maxLabelSize)
        self.connect(self.slidervthB2, SIGNAL('sliderReleased()'), self.onSlidervthB2)
        self.connect(self.sliderValuevthB2, SIGNAL('clicked()'), self.onSliderValuevthB2)

        sliderLabel = QLabel('v [km/s]')
        self.slidervB2 = FloatSlider(Qt.Horizontal)
        self.sliderLeftvB2 = ExtendedQLabel('0')
        self.sliderRightvB2 = ExtendedQLabel('180')
        self.slidervB2.setMinimum(-10)
        self.slidervB2.setMaximum(10)
        self.slidervB2.setValue(self.dopplerVelocity2Input)
        self.slidervB2.setTracking(True)
        self.slidervB2.setTickPosition(QSlider.TicksBothSides)
        self.slidervB2.setMinimumWidth(self.minSizeSlider)
        self.slidervB2.setMaximumWidth(self.maxSizeSlider)
        self.sliderValuevB2 = ExtendedQLabel('{0:6.1f}'.format(self.dopplerVelocity2Input))
        for l, w in enumerate([sliderLabel, self.sliderLeftvB2, self.slidervB2, self.sliderRightvB2, self.sliderValuevB2]):
            boxB2.addWidget(w, 4, l)
            boxB2.setAlignment(w, Qt.AlignVCenter)
        self.sliderValuevB2.setMinimumWidth(self.minLabelSize)
        self.sliderValuevB2.setMaximumWidth(self.maxLabelSize)
        self.connect(self.slidervB2, SIGNAL('sliderReleased()'), self.onSlidervB2)
        self.connect(self.sliderValuevB2, SIGNAL('clicked()'), self.onSliderValuevB2)

        sliderLabel = QLabel('tau')
        self.slidertau2 = FloatSlider(Qt.Horizontal)
        self.sliderLefttau2 = ExtendedQLabel('0')
        self.sliderRighttau2 = ExtendedQLabel('180')
        self.slidertau2.setMinimum(0.0)
        self.slidertau2.setMaximum(10.0)
        self.slidertau2.setValue(self.tau2Input)
        self.slidertau2.setTracking(True)
        self.slidertau2.setTickPosition(QSlider.TicksBothSides)
        self.slidertau2.setMinimumWidth(self.minSizeSlider)
        self.slidertau2.setMaximumWidth(self.maxSizeSlider)
        self.sliderValuetau2 = ExtendedQLabel('{0:6.1f}'.format(self.tau2Input))
        for l, w in enumerate([sliderLabel, self.sliderLefttau2, self.slidertau2, self.sliderRighttau2, self.sliderValuetau2]):
            boxB2.addWidget(w, 5, l)
            boxB2.setAlignment(w, Qt.AlignVCenter)
        self.sliderValuetau2.setMinimumWidth(self.minLabelSize)
        self.sliderValuetau2.setMaximumWidth(self.maxLabelSize)
        self.connect(self.slidertau2, SIGNAL('sliderReleased()'), self.onSlidertau2)
        self.connect(self.sliderValuetau2, SIGNAL('clicked()'), self.onSliderValuetau2)


        comp1Group = QGroupBox("Component 1")
        comp1Group.setLayout(boxB1)
        vboxB1 = QVBoxLayout()
        vboxB1.addWidget(comp1Group)

        comp2Group = QGroupBox("Component 2")
        comp2Group.setLayout(boxB2)
        vboxB2 = QVBoxLayout()
        vboxB2.addWidget(comp2Group)


        vboxB = QVBoxLayout()
        vboxB.addLayout(vboxB1)
        vboxB.addLayout(vboxB2)

        # Radiative transfer
        # 
        radTranGroup = QGroupBox("Radiative transfer")
        self.radioThin = QRadioButton("Optically thin")
        self.radioExact = QRadioButton("Exact")
        self.radioExact.setChecked(True)
        self.radioAtompol = QRadioButton("Exact Without Atompol")

        normalizationGroup = QGroupBox("Normalization")
        self.radioMaxNorm = QRadioButton("Normalize maximum")
        self.radioPeakNorm = QRadioButton("Normalize peak")
        self.radioMaxNorm.setChecked(True)

        self.connect(self.radioThin, SIGNAL('clicked()'), self.onRadioThin)
        self.connect(self.radioExact, SIGNAL('clicked()'), self.onRadioExact)
        self.connect(self.radioAtompol, SIGNAL('clicked()'), self.onRadioAtompol)

        self.connect(self.radioMaxNorm, SIGNAL('clicked()'), self.onRadioMaxNorm)
        self.connect(self.radioPeakNorm, SIGNAL('clicked()'), self.onRadioPeakNorm)

        
        vboxRadTran = QHBoxLayout()
        vboxRadTran.addWidget(self.radioThin)
        vboxRadTran.addWidget(self.radioExact)
        vboxRadTran.addWidget(self.radioAtompol)        
        vboxRadTran.addStretch(1)
        radTranGroup.setLayout(vboxRadTran)        

        vboxRadTran = QHBoxLayout()
        vboxRadTran.addWidget(self.radioMaxNorm)
        vboxRadTran.addWidget(self.radioPeakNorm)
        vboxRadTran.addStretch(1)
        normalizationGroup.setLayout(vboxRadTran)

        boundary = QGridLayout()
        self.checkAllen = QCheckBox("Use Allen")
        self.checkAllen.toggle()
        self.connect(self.checkAllen, SIGNAL('stateChanged(int)'), self.onCheckAllen)        

        boundary.addWidget(self.checkAllen, 0, 0)

        title = QLabel('I0:')
        self.I0 = QLineEdit("0")
        boundary.addWidget(title, 0, 1)
        boundary.addWidget(self.I0, 0, 2)
        self.connect(self.I0, SIGNAL('editingFinished()'), self.onChangeI0)        

        title = QLabel('Q0:')
        self.Q0 = QLineEdit("0")        
        boundary.addWidget(title, 0, 3)
        boundary.addWidget(self.Q0, 0, 4)
        self.connect(self.Q0, SIGNAL('editingFinished()'), self.onChangeQ0)

        title = QLabel('U0:')
        self.U0 = QLineEdit("0")
        boundary.addWidget(title, 0, 5)
        boundary.addWidget(self.U0, 0, 6)
        self.connect(self.U0, SIGNAL('editingFinished()'), self.onChangeU0)

        title = QLabel('V0:')
        self.V0 = QLineEdit("0")
        boundary.addWidget(title, 0, 7)
        boundary.addWidget(self.V0, 0, 8)
        self.connect(self.V0, SIGNAL('editingFinished()'), self.onChangeV0)

        self.I0.setEnabled(False)
        self.Q0.setEnabled(False)
        self.U0.setEnabled(False)
        self.V0.setEnabled(False)

        title = QLabel('N slabs:')
        self.nslabs = QComboBox()
        self.nslabs.addItem("One")
        self.nslabs.addItem("Two vertical")
        self.nslabs.addItem("Two horizontal")
        boundary.addWidget(title, 0, 9)
        boundary.addWidget(self.nslabs, 0, 10)
        self.connect(self.nslabs, SIGNAL('activated(int)'), self.onChangenSlabs)
        

        wave = QGridLayout()
        # title = QLabel('Allen:')
        # self.allen = QLineEdit("0")
        # wave.addWidget(title, 0, 0)
        # wave.addWidget(self.allen, 0, 1)

        title = QLabel('wl:')
        self.leftWave = QLineEdit(str(self.lambdaAxisInput[0]))
        wave.addWidget(title, 0, 2)
        wave.addWidget(self.leftWave, 0, 3)
        self.connect(self.leftWave, SIGNAL('editingFinished()'), self.onChangeleftWave)

        title = QLabel('wr:')
        self.rightWave = QLineEdit(str(self.lambdaAxisInput[1]))
        wave.addWidget(title, 0, 4)
        wave.addWidget(self.rightWave, 0, 5)
        self.connect(self.rightWave, SIGNAL('editingFinished()'), self.onChangerightWave)

        title = QLabel('step:')
        self.stepWave = QLineEdit(str(self.nLambdaInput))
        wave.addWidget(title, 0, 6)
        wave.addWidget(self.stepWave, 0, 7)
        self.connect(self.stepWave, SIGNAL('editingFinished()'), self.onChangestepWave)

        vboxR = QVBoxLayout()        
        vboxR.addWidget(radTranGroup)
        vboxR.addWidget(normalizationGroup)
        vboxR.addLayout(boundary)
        vboxR.addLayout(wave)


        # Radiative transfer
        # 
        boxO1 = QGridLayout()
        sliderLabel = QLabel('theta [deg]')
        self.slidertheta = FloatSlider(Qt.Horizontal)
        self.slidertheta.setMinimum(0.0)
        self.slidertheta.setMaximum(180)
        self.slidertheta.setValue(self.anglesInput[0])
        self.slidertheta.setTracking(True)
        self.slidertheta.setTickPosition(QSlider.TicksBothSides)
        self.slidertheta.setMinimumWidth(self.minSizeSlider)
        self.slidertheta.setMaximumWidth(self.maxSizeSlider)
        self.sliderValuetheta = ExtendedQLabel('{0:6.1f}'.format(self.anglesInput[0]))
        for l, w in enumerate([sliderLabel, self.slidertheta, self.sliderValuetheta]):
            boxO1.addWidget(w, 0, l)
            boxO1.setAlignment(w, Qt.AlignVCenter)
        self.sliderValuetheta.setMinimumWidth(self.minLabelSize)
        self.sliderValuetheta.setMaximumWidth(self.maxLabelSize)
        self.connect(self.slidertheta, SIGNAL('sliderReleased()'), self.onSliderTheta)
        self.connect(self.sliderValuetheta, SIGNAL('clicked()'), self.onSliderValueTheta)
        if (self.checkAllen.isChecked()):
            theta = float(self.sliderValuetheta.text())
            mu = np.cos(theta * np.pi / 180.0)
            i0 = i0Allen.i0Allen(10830.0, mu)
            self.I0.setText(str(i0))
            self.boundaryInput[0] = i0


        sliderLabel = QLabel('phi [deg]')
        self.sliderphi = FloatSlider(Qt.Horizontal)
        self.sliderphi.setMinimum(0.0)
        self.sliderphi.setMaximum(180)
        self.sliderphi.setValue(self.anglesInput[1])
        self.sliderphi.setTracking(True)
        self.sliderphi.setTickPosition(QSlider.TicksBothSides)
        self.sliderphi.setMinimumWidth(self.minSizeSlider)
        self.sliderphi.setMaximumWidth(self.maxSizeSlider)
        self.sliderValuephi = ExtendedQLabel('{0:6.1f}'.format(self.anglesInput[1]))
        for l, w in enumerate([sliderLabel, self.sliderphi, self.sliderValuephi]):
            boxO1.addWidget(w, 1, l)
            boxO1.setAlignment(w, Qt.AlignVCenter)
        self.sliderValuephi.setMinimumWidth(self.minLabelSize)
        self.sliderValuephi.setMaximumWidth(self.maxLabelSize)
        self.connect(self.sliderphi, SIGNAL('sliderReleased()'), self.onSliderPhi)
        self.connect(self.sliderValuephi, SIGNAL('clicked()'), self.onSliderValuePhi)

        sliderLabel = QLabel('gamma [deg]')
        self.slidergamma = FloatSlider(Qt.Horizontal)
        self.slidergamma.setMinimum(0)
        self.slidergamma.setMaximum(180)
        self.slidergamma.setValue(self.anglesInput[2])
        self.slidergamma.setTracking(True)
        self.slidergamma.setTickPosition(QSlider.TicksBothSides)
        self.slidergamma.setMinimumWidth(self.minSizeSlider)
        self.slidergamma.setMaximumWidth(self.maxSizeSlider)
        self.sliderValuegamma = ExtendedQLabel('{0:6.1f}'.format(self.anglesInput[2]))
        for l, w in enumerate([sliderLabel, self.slidergamma, self.sliderValuegamma]):
            boxO1.addWidget(w, 2, l)
            boxO1.setAlignment(w, Qt.AlignVCenter)
        self.sliderValuegamma.setMinimumWidth(self.minLabelSize)
        self.sliderValuegamma.setMaximumWidth(self.maxLabelSize)
        self.connect(self.slidergamma, SIGNAL('sliderReleased()'), self.onSliderGamma)
        self.connect(self.sliderValuegamma, SIGNAL('clicked()'), self.onSliderValueGamma)

        sliderLabel = QLabel('height [deg]')
        self.sliderheight = FloatSlider(Qt.Horizontal)
        self.sliderheight.setMinimum(0.1)
        self.sliderheight.setMaximum(20.0)
        self.sliderheight.setValue(self.hInput)
        self.sliderheight.setTracking(True)
        self.sliderheight.setTickPosition(QSlider.TicksBothSides)
        self.sliderheight.setMinimumWidth(self.minSizeSlider)
        self.sliderheight.setMaximumWidth(self.maxSizeSlider)
        self.sliderValueheight = ExtendedQLabel('{0:6.1f}'.format(self.hInput))
        for l, w in enumerate([sliderLabel, self.sliderheight, self.sliderValueheight]):
            boxO1.addWidget(w, 3, l)
            boxO1.setAlignment(w, Qt.AlignVCenter)
        self.sliderValueheight.setMinimumWidth(self.minLabelSize)
        self.sliderValueheight.setMaximumWidth(self.maxLabelSize)
        self.connect(self.sliderheight, SIGNAL('sliderReleased()'), self.onSliderHeight)
        self.connect(self.sliderValueheight, SIGNAL('clicked()'), self.onSliderValueHeight)

        sliderLabel = QLabel('a')
        self.sliderDamp = FloatSlider(Qt.Horizontal)
        self.sliderDamp.setMinimum(0.0)
        self.sliderDamp.setMaximum(2.0)
        self.sliderDamp.setValue(self.dampingInput)
        self.sliderDamp.setTracking(True)
        self.sliderDamp.setTickPosition(QSlider.TicksBothSides)
        self.sliderDamp.setMinimumWidth(self.minSizeSlider)
        self.sliderDamp.setMaximumWidth(self.maxSizeSlider)
        self.sliderValueDamp = ExtendedQLabel('{0:6.1f}'.format(self.dampingInput))
        for l, w in enumerate([sliderLabel, self.sliderDamp, self.sliderValueDamp]):
            boxO1.addWidget(w, 4, l)
            boxO1.setAlignment(w, Qt.AlignVCenter)
        self.sliderValueDamp.setMinimumWidth(self.minLabelSize)
        self.sliderValueDamp.setMaximumWidth(self.maxLabelSize)
        self.connect(self.sliderDamp, SIGNAL('sliderReleased()'), self.onSliderDamp)
        self.connect(self.sliderValueDamp, SIGNAL('clicked()'), self.onSliderValueDamp)

        sliderLabel = QLabel('ff')
        self.sliderff = FloatSlider(Qt.Horizontal)
        self.sliderff.setMinimum(0.0)
        self.sliderff.setMaximum(1.0)
        self.sliderff.setValue(self.ffInput)
        self.sliderff.setTracking(True)
        self.sliderff.setTickPosition(QSlider.TicksBothSides)
        self.sliderff.setMinimumWidth(self.minSizeSlider)
        self.sliderff.setMaximumWidth(self.maxSizeSlider)
        self.sliderValueff = ExtendedQLabel('{0:6.1f}'.format(self.ffInput))
        for l, w in enumerate([sliderLabel, self.sliderff, self.sliderValueff]):
            boxO1.addWidget(w, 5, l)
            boxO1.setAlignment(w, Qt.AlignVCenter)
        self.sliderValueff.setMinimumWidth(self.minLabelSize)
        self.sliderValueff.setMaximumWidth(self.maxLabelSize)
        self.connect(self.sliderff, SIGNAL('sliderReleased()'), self.onSliderff)
        self.connect(self.sliderValueff, SIGNAL('clicked()'), self.onSliderValueff)

        sliderLabel = QLabel('S2/S1 s1')
        self.sliderS2S1 = FloatSlider(Qt.Horizontal)
        self.sliderS2S1.setMinimum(0.0)
        self.sliderS2S1.setMaximum(2.0)
        self.sliderS2S1.setValue(self.betaInput)
        self.sliderS2S1.setTracking(True)
        self.sliderS2S1.setTickPosition(QSlider.TicksBothSides)
        self.sliderS2S1.setMinimumWidth(self.minSizeSlider)
        self.sliderS2S1.setMaximumWidth(self.maxSizeSlider)
        self.sliderValueS2S1 = ExtendedQLabel('{0:6.1f}'.format(self.betaInput))
        for l, w in enumerate([sliderLabel, self.sliderS2S1, self.sliderValueS2S1]):
            boxO1.addWidget(w, 6, l)
            boxO1.setAlignment(w, Qt.AlignVCenter)
        self.sliderValueS2S1.setMinimumWidth(self.minLabelSize)
        self.sliderValueS2S1.setMaximumWidth(self.maxLabelSize)
        self.connect(self.sliderS2S1, SIGNAL('sliderReleased()'), self.onSliderS2S1)
        self.connect(self.sliderValueS2S1, SIGNAL('clicked()'), self.onSliderValueS2S1)

        sliderLabel = QLabel('S2/S1 s2')
        self.sliderS2S1B = FloatSlider(Qt.Horizontal)
        self.sliderS2S1B.setMinimum(0.0)
        self.sliderS2S1B.setMaximum(2.0)
        self.sliderS2S1B.setValue(self.betaInput)
        self.sliderS2S1B.setTracking(True)
        self.sliderS2S1B.setTickPosition(QSlider.TicksBothSides)
        self.sliderS2S1B.setMinimumWidth(self.minSizeSlider)
        self.sliderS2S1B.setMaximumWidth(self.maxSizeSlider)
        self.sliderValueS2S1B = ExtendedQLabel('{0:6.1f}'.format(self.beta2Input))
        for l, w in enumerate([sliderLabel, self.sliderS2S1B, self.sliderValueS2S1B]):
            boxO1.addWidget(w, 7, l)
            boxO1.setAlignment(w, Qt.AlignVCenter)
        self.sliderValueS2S1B.setMinimumWidth(self.minLabelSize)
        self.sliderValueS2S1B.setMaximumWidth(self.maxLabelSize)
        self.connect(self.sliderS2S1B, SIGNAL('sliderReleased()'), self.onSliderS2S1B)
        self.connect(self.sliderValueS2S1B, SIGNAL('clicked()'), self.onSliderValueS2S1B)

        sliderLabel = QLabel('multiplet')
        self.slidermultiplet = QSlider(Qt.Horizontal)
        self.slidermultiplet.setRange(1, 4)
        self.slidermultiplet.setValue(self.transInput)
        self.slidermultiplet.setTracking(True)
        self.slidermultiplet.setTickPosition(QSlider.TicksBothSides)
        self.slidermultiplet.setMinimumWidth(self.minSizeSlider)
        self.slidermultiplet.setMaximumWidth(self.maxSizeSlider)
        self.sliderValuemultiplet = ExtendedQLabel(self.multiplets[self.transInput-1]) 
        for l, w in enumerate([sliderLabel, self.slidermultiplet, self.sliderValuemultiplet]):
            boxO1.addWidget(w, 8, l)
            boxO1.setAlignment(w, Qt.AlignVCenter)
        self.sliderValuemultiplet.setMinimumWidth(self.minLabelSize)
        self.sliderValuemultiplet.setMaximumWidth(self.maxLabelSize)
        self.connect(self.slidermultiplet, SIGNAL('sliderReleased()'), self.onSlidermultiplet)


        obsGroup = QGroupBox("Observation")
        obsGroup.setLayout(boxO1)
        vboxO1 = QVBoxLayout()
        vboxO1.addWidget(obsGroup)
        # vboxO1.addWidget(self.allenLabel)

        
        self.calculateButton = QPushButton("Calculate")
        self.loadButton = QPushButton("Load observation")
        self.resetButton = QPushButton("Reset observation")
        if (self.obsFile != ''):
            self.loadedFile = ExtendedQLabel('Loaded: {0}'.format(self.obsFile))
        else:
            self.loadedFile = ExtendedQLabel('Loaded: none')

        boxO2 = QVBoxLayout()
        boxO2.addWidget(self.calculateButton)
        boxO2.addWidget(self.loadButton)
        boxO2.addWidget(self.resetButton)
        boxO2.addWidget(self.loadedFile)

        calcGroup = QGroupBox("Calculate")
        calcGroup.setLayout(boxO2)
        vboxO2 = QVBoxLayout()
        vboxO2.addWidget(calcGroup)      

        vboxO = QVBoxLayout()
        vboxO.addLayout(vboxO1)
        vboxO.addLayout(vboxO2)

        self.connect(self.calculateButton, SIGNAL('clicked()'), self.redrawProfiles)
        self.connect(self.loadButton, SIGNAL('clicked()'), self.loadObservation)
        self.connect(self.resetButton, SIGNAL('clicked()'), self.resetObservation)
        
        # Final layout
        #
        hbox1 = QHBoxLayout()
        hbox1.addLayout(layout)
        hbox1.addLayout(vboxB)
        hbox1.addLayout(vboxO)

        hbox2 = QHBoxLayout()
        hbox2.addLayout(vboxR)        

        vbox = QVBoxLayout()
        vbox.addLayout(hbox1)
        vbox.addLayout(hbox2)

        if (self.obsFile != ''):
            self.obsStokes = np.loadtxt(self.obsFile)


        # self.textbox = QLineEdit()
        # self.textbox.setMinimumWidth(200)
        # self.connect(self.textbox, SIGNAL('editingFinished ()'), self.on_draw)
        
        # self.draw_button = QPushButton("&Draw")
        # self.connect(self.draw_button, SIGNAL('clicked()'), self.on_draw)
        
        # self.grid_cb = QCheckBox("Show &Grid")
        # self.grid_cb.setChecked(False)
        # self.connect(self.grid_cb, SIGNAL('stateChanged(int)'), self.on_draw)
        
        # slider_label = QLabel('Bar width (%):')
        # self.slider = QSlider(Qt.Horizontal)
        # self.slider.setRange(1, 100)
        # self.slider.setValue(20)
        # self.slider.setTracking(True)
        # self.slider.setTickPosition(QSlider.TicksBothSides)
        # self.connect(self.slider, SIGNAL('valueChanged(int)'), self.on_draw)
        
        # #
        # # Layout with box sizers
        # # 
        # hbox = QHBoxLayout()
        
        # for w in [  self.textbox, self.draw_button, self.grid_cb,
        #             slider_label, self.slider]:
        #     hbox.addWidget(w)
        #     hbox.setAlignment(w, Qt.AlignVCenter)
        
        # vbox = QVBoxLayout()
        # vbox.addWidget(self.canvas)
        # vbox.addWidget(self.mpl_toolbar)
        # vbox.addLayout(hbox)
        # vbox.addLayout(vboxB1)
        
        self.main_frame.setLayout(vbox)
        self.setCentralWidget(self.main_frame)
    
    def create_status_bar(self):
        self.status_text = QLabel("OK")
        self.statusBar().addWidget(self.status_text, 1)

    def setMultitermHe(self):
        pass

    def setMultitermNa(self):
        pass

    def on_increasefont(self):
        self.fontSize += 1
        self.font = QFont("SansSerif", self.fontSize)

        QApplication.setFont(self.font)

    def on_decreasefont(self):
        self.fontSize -= 1
        self.font = QFont("SansSerif", self.fontSize)

        QApplication.setFont(self.font)
        
    def create_menu(self):        
        """Create main menu
        
        Returns:
            TYPE: None
        """
        self.file_menu = self.menuBar().addMenu("&Multiterm")

        he_action = self.create_action("&He", slot=self.setMultitermHe, 
            shortcut="Ctrl+H", tip="Change to He I atom model")

        na_action = self.create_action("&Na", slot=self.setMultitermNa, 
            shortcut="Ctrl+N", tip="Change to Na I atom model")
        
        quit_action = self.create_action("&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")
        
        self.add_actions(self.file_menu, (he_action, na_action, None, quit_action))

        self.font_menu = self.menuBar().addMenu("&Font")
        increasefont_action = self.create_action("&Larger", 
            shortcut='Ctrl++', slot=self.on_increasefont, 
            tip='Increase font size')

        decreasefont_action = self.create_action("&Smaller", 
            shortcut='Ctrl+-', slot=self.on_decreasefont, 
            tip='Decrease font size')

        self.add_actions(self.font_menu, (increasefont_action, decreasefont_action))
        
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About", 
            shortcut='F1', slot=self.on_about, 
            tip='About hazel')
        
        self.add_actions(self.help_menu, (about_action,))

    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(  self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False, 
                        signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action

    def closeEvent(self, event):
        print("Saving state...")
        self.saveConfig()


def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()
