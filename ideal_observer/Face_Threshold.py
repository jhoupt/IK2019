#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v3.0.4),
    on Thu Mar 14 21:00:03 2019
If you publish work using this script please cite the PsychoPy publications:
    Peirce, JW (2007) PsychoPy - Psychophysics software in Python.
        Journal of Neuroscience Methods, 162(1-2), 8-13.
    Peirce, JW (2009) Generating stimuli for neuroscience using PsychoPy.
        Frontiers in Neuroinformatics, 2:10. doi: 10.3389/neuro.11.010.2008
"""

from __future__ import absolute_import, division
from psychopy import locale_setup, sound, gui, visual, core, data, event, logging, clock
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)
import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle
import os  # handy system and path functions
import sys  # to get file system encoding
from PIL import Image
from scipy.stats import norm
import csv

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)

# Store info about the experiment session
psychopyVersion = '3.0.4'
expName = 'Face-Threshold'  # from the Builder filename that created this script
expInfo = {'Face': ['Real', 'Schematic'], 'participant': ''}
#expInfo.addField('Face:', choices=['Real', 'Schematic'])
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName)
if dlg.OK == False:
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
expInfo['psychopyVersion'] = psychopyVersion

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath='/Users/joe/Box/Teaching/IK Mathematical Cognitive Modeling/Examples/Face-Threshold.py',
    savePickle=True, saveWideText=True,
    dataFileName=filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp

# Start Code - component code to be run before the window creation

# Setup the Window
win = visual.Window(
    size=(1024, 768), fullscr=True, screen=0,
    allowGUI=False, allowStencil=False,
    monitor='testMonitor', color=[0,0,0], colorSpace='rgb',
    blendMode='avg', useFBO=True,
    units='height')
# store frame rate of monitor if we can measure it
expInfo['frameRate'] = win.getActualFrameRate()
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess

# Initialize components for Routine "trial"
trialClock = core.Clock()
fixation_cross = visual.ShapeStim(
    win=win, name='fixation_cross', vertices='cross',
    size=(0.01, 0.01),
    ori=0, pos=(0, 0),
    lineWidth=.001, lineColor=[1,1,1], lineColorSpace='rgb',
    fillColor=[1,1,1], fillColorSpace='rgb',
    opacity=1, depth=-1.0, interpolate=True)

# Load Stimuli
if expInfo['Face'] == 'Real' :
    happy_face = Image.open("Happy_real.bmp")
    sad_face = Image.open("Sad_real.bmp")
else :
    happy_face = Image.open("Happy_schematic.bmp")
    sad_face = Image.open("Sad_schematic.bmp")
happy_array = np.array(happy_face) / 255 - 0.5
sad_array = np.array(sad_face) / 255 - 0.5

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

# --------Prepare to start Staircase "trials" --------
# set up handler to look after next chosen value etc
if expInfo['participant'] == 'ideal' :
    trials = data.QuestHandler(startVal = 0.5, 
                      startValSd = 0.2, delta=0,
                      pThreshold = 0.63, 
                      nTrials=1E3, minVal=0, maxVal=1)
else :
    trials = data.QuestHandler(startVal = 0.5, 
                      startValSd = 0.2, 
                      pThreshold = 0.63, 
                      nTrials=50, minVal=0, maxVal=1)
   
thisExp.addLoop(trials)  # add the loop to the experiment
level = thisTrial = 0.5  # initialise some vals

for thisTrial in trials:
    currentLoop = trials
    level = thisTrial
    # update component parameters for each repeat
    key_resp_2 = event.BuilderKeyResponse()
    
    if np.random.random() < .5 :
        trial_array = happy_array * level + np.random.normal(0, .2, happy_array.shape)
        correct_key = str('q')
    else :
        trial_array = sad_array * level + np.random.normal(0, .2, sad_array.shape)
        correct_key = str('p')

    trial_array[trial_array > .5] = .5; trial_array[trial_array < -.5] = -.5

    if  expInfo['participant'] == 'ideal' :
        log_likelihood_happy = np.sum(norm.logpdf(trial_array - level*happy_array, loc=0, scale=.2))
        log_likelihood_sad = np.sum(norm.logpdf(trial_array - level*sad_array, loc=0, scale=.2))
        choose_happy = log_likelihood_happy > log_likelihood_sad
        if (choose_happy and correct_key == 'q') or  (not choose_happy and correct_key == 'p'):
           key_resp_2.corr = 1
        else :
           key_resp_2.corr = 0
        key_resp_2.rt = 0
    else :
        trial_array = Image.fromarray((trial_array + 0.5) * 255)
        face = visual.ImageStim(win, trial_array)
        # ------Prepare to start Routine "trial"-------
        t = 0
        trialClock.reset()  # clock
        frameN = -1
        continueRoutine = True
        routineTimer.add(2.000000)
        
        # keep track of which components have finished
        trialComponents = [face, fixation_cross, key_resp_2]
        for thisComponent in trialComponents:
            if hasattr(thisComponent, 'status'):
                thisComponent.status = NOT_STARTED
        
        # -------Start Routine "trial"-------
        while continueRoutine and routineTimer.getTime() > 0:
            # get current time
            t = trialClock.getTime()
            frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
            # update/draw components on each frame
            
            # *face* updates
            if t >= 0.4 and face.status == NOT_STARTED:
                # keep track of start time/frame for later
                face.tStart = t
                face.frameNStart = frameN  # exact frame index
                face.setAutoDraw(True)
            frameRemains = 0.4 + .25- win.monitorFramePeriod * 0.75  # most of one frame period left
            if face.status == STARTED and t >= frameRemains:
                face.setAutoDraw(False)
            
            # *fixation_cross* updates
            if t >= 0.0 and fixation_cross.status == NOT_STARTED:
                # keep track of start time/frame for later
                fixation_cross.tStart = t
                fixation_cross.frameNStart = frameN  # exact frame index
                fixation_cross.setAutoDraw(True)
            frameRemains = 0.0 + .15- win.monitorFramePeriod * 0.75  # most of one frame period left
            if fixation_cross.status == STARTED and t >= frameRemains:
                fixation_cross.setAutoDraw(False)
            
            # *key_resp_2* updates
            if t >= 0.0 and key_resp_2.status == NOT_STARTED:
                # keep track of start time/frame for later
                key_resp_2.tStart = t
                key_resp_2.frameNStart = frameN  # exact frame index
                key_resp_2.status = STARTED
                # keyboard checking is just starting
                win.callOnFlip(key_resp_2.clock.reset)  # t=0 on next screen flip
                event.clearEvents(eventType='keyboard')
            frameRemains = 0.0 + 3- win.monitorFramePeriod * 0.75  # most of one frame period left
            if key_resp_2.status == STARTED and t >= frameRemains:
                key_resp_2.status = FINISHED
            if key_resp_2.status == STARTED:
                theseKeys = event.getKeys(keyList=['q', 'p', 'esc'])
                
                # check for quit:
                if "escape" in theseKeys:
                    endExpNow = True
                if len(theseKeys) > 0:  # at least one key was pressed
                    key_resp_2.keys = theseKeys[-1]  # just the last key pressed
                    key_resp_2.rt = key_resp_2.clock.getTime()
                    # was this 'correct'?
                    if (key_resp_2.keys == correct_key):
                        key_resp_2.corr = 1
                    else:
                        key_resp_2.corr = 0
                    # a response ends the routine
                    continueRoutine = False
            
            # check for quit (typically the Esc key)
            if endExpNow or event.getKeys(keyList=["escape"]):
                core.quit()
            
            # check if all components have finished
            if not continueRoutine:  # a component has requested a forced-end of Routine
                break
            continueRoutine = False  # will revert to True if at least one component still running
            for thisComponent in trialComponents:
                if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                    continueRoutine = True
                    break  # at least one component has not yet finished
            
            # refresh the screen
            if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
                win.flip()
        
        # -------Ending Routine "trial"-------
        if t < 2 :
            win.flip()
            core.wait(2 - t)
            
        for thisComponent in trialComponents:
            if hasattr(thisComponent, "setAutoDraw"):
                thisComponent.setAutoDraw(False)

    # store data for trials (StairHandler)
    trials.addResponse(key_resp_2.corr)
    trials.addOtherData('key_resp_2.rt', key_resp_2.rt)
    thisExp.nextEntry()
    
# staircase completed
print('The threshold is {s}'.format(s=trials.mean()))
with open(filename + '_threshold.csv', 'w') as f:    f.write("Threshold,Face,participant,date,expName\n")
    f.write(str(trials.mean()) + ',' + expInfo['Face']  + ',' + expInfo['participant'] + ',' + expInfo['date'] + ',' + expName)


# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename+'.csv')
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
