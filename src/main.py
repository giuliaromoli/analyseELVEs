# created:15/02/2023

import yaml
import ROOT
import os
from easydict import EasyDict as ed
from src.utils.getWindowData import getWindowLimits
from src.utils.getWindowData import getData

ROOT.gROOT.SetBatch(ROOT.kTRUE)  # do not open figures!

if __name__ == '__main__':

    listEvents = []
    for inputFile in os.listdir(r'/home/jule/Scrivania/Mini-Euso/analyseELVEs/src/config'):
        if inputFile.endswith('.yaml'):
            listEvents.append(inputFile.split('.')[0].split('config-')[1])
    #print(listEvents)
    for ee in [0]:#range(len(listEvents)):
        #eventLabel = listEvents[ee]
        eventLabel = "20-01-2021"
        eventPath = r'/home/jule/Scrivania/Mini-Euso/analyseELVEs/results/' + eventLabel + '/'

        if not os.path.exists(eventPath):
            os.makedirs(eventPath)

        print('start analysing ' + eventLabel)

        with open(r'/home/jule/Scrivania/Mini-Euso/analyseELVEs/src/config/config-' + eventLabel + '.yaml') as file:
            config = yaml.safe_load(file)
        configFile = ed(config)

        # get window limits:
        if not os.path.exists(eventPath + "/infoWindow/infoWindow.txt"):
            (startWindow,endWindow) = getWindowLimits(eventPath, configFile)
        else:
            with open(eventPath + "/infoWindow/infoWindow.txt") as f:
                lines = f.readlines()
                startWindow = int(lines[0].split(':')[1])
                endWindow = int(lines[1].split(':')[1])
        #print(startWindow,endWindow)
        #get data:
        photonCounterDataWindow = getData(startWindow,endWindow,eventPath,configFile)