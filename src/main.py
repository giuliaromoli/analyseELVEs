# created:15/02/2023

import yaml
import ROOT
import os
from easydict import EasyDict as ed
from src.utils.getWindowData import getWindowLimits
from src.utils.getWindowData import getData
from src.utils.getPixelsInfos import getPixelsInfos
from src.utils.getPixelsInfos import getPhiEnergy
from src.utils.getPixelsInfos import getPixelsPlots
from src.utils.getCircleCenter import getCircleCenter
from src.utils.fillHisto import getPolarHisto
from src.utils.getPixelsInfosII import getPixelsInfosII


ROOT.gROOT.SetBatch(ROOT.kTRUE)  # do not open figures!

if __name__ == '__main__':

    listEvents = []
    for inputFile in os.listdir(r'/home/jule/Scrivania/Mini-Euso/analyseELVEs/src/config'):
        if inputFile.endswith('.yaml'):
            listEvents.append(inputFile.split('.')[0].split('config-')[1])
    # print(listEvents)
    for ee in [0]:  # range(len(listEvents)):
        eventLabel = "09-01-2021n1"  # listEvents[ee]
        eventPath = r'/home/jule/Scrivania/Mini-Euso/analyseELVEs/results/' + eventLabel + '/'

        if not os.path.exists(eventPath):
            os.makedirs(eventPath)

        print('start analysing ' + eventLabel)

        with open(r'/home/jule/Scrivania/Mini-Euso/analyseELVEs/src/config/config-' + eventLabel + '.yaml') as file:
            config = yaml.safe_load(file)
        configFile = ed(config)

        # get window limits:
        if not os.path.exists(eventPath + "/infoWindow/infoWindow.txt"):
            (startWindow, endWindow) = getWindowLimits(eventPath, configFile)
        else:
            with open(eventPath + "/infoWindow/infoWindow.txt") as f:
                lines = f.readlines()
                startWindow = int(lines[0].split(':')[1])
                endWindow = int(lines[1].split(':')[1])
        #print(startWindow,endWindow)

        # get data:
        photonCounterDataWindow = getData(startWindow, endWindow, eventPath, configFile)
        #continue
        # fit pixels with GaussExp:
        if not os.path.exists(eventPath + "/pixels/infoPixels.txt"):
            pixelsInfos = getPixelsInfos(eventPath + "/pixels", configFile, photonCounterDataWindow, startWindow)
        else:
            pixelsInfos = []
            with open(eventPath + "/pixels/infoPixels.txt") as f:
                for line in f:
                    splittedLine = line.split('\t')
                    pixelLabel = int(splittedLine[0])
                    pixX = int(splittedLine[1])
                    pixY = int(splittedLine[2])
                    x = float(splittedLine[3])
                    y = float(splittedLine[4])
                    nPeak = int(splittedLine[5])
                    timeMax = int(splittedLine[6])
                    maxEnergy = float(splittedLine[7])
                    sigmaMaxEnergy = float(splittedLine[8])
                    fullWidth = int(splittedLine[9])
                    deExcTime = int(splittedLine[10])
                    totEnergy = float(splittedLine[11])
                    sigmaTotEnergy = float(splittedLine[12])
                    energyPeak = float(splittedLine[13])
                    sigmaEnergyPeak = float(splittedLine[14])
                    # from fit with gaussExp:
                    tauFit = float(splittedLine[15])
                    sigmaFit = float(splittedLine[16])
                    meanFit = int(splittedLine[17])
                    areaFit = float(splittedLine[18])
                    maxEnergyFit = float(splittedLine[19])
                    sigmaMaxEnergyFit = float(splittedLine[20])
                    energyPeakFit = float(splittedLine[21])
                    sigmaEnergyPeakFit = float(splittedLine[22])
                    thresholdSinglePeak = float(splittedLine[23])
                    lambdaFit = float(splittedLine[24])
                    totDuration = int(splittedLine[25])
                    totDeExcTime = int(splittedLine[26].split('\n')[0])
                    pixelsInfos.append((pixelLabel, pixX, pixY, x, y, nPeak, timeMax, maxEnergy, sigmaMaxEnergy,
                                        fullWidth, deExcTime, totEnergy, sigmaTotEnergy, energyPeak, sigmaEnergyPeak,
                                        tauFit, sigmaFit, meanFit, areaFit, maxEnergyFit, sigmaMaxEnergyFit,
                                        energyPeakFit, sigmaEnergyPeakFit, thresholdSinglePeak, lambdaFit, totDuration,
                                        totDeExcTime))

        #getPixelsPlots(eventPath + "/pixelsPlots", pixelsInfos, configFile)
        #continue
        #fit circle:
        if not os.path.exists(eventPath + "/circleFrames/infoCenter.txt"):
            (bestXCircleCenter, bestYCircleCenter, errCentreX, errCentreY) = getCircleCenter(eventPath + "/circleFrames/", pixelsInfos, photonCounterDataWindow, configFile, startWindow)
        else:
            with open(eventPath + "/circleFrames/infoCenter.txt") as f:
                lines = f.readlines()
                bestXCircleCenter = float(lines[1].split(':')[1].split(' ')[1])
                bestYCircleCenter = float(lines[3].split(':')[1].split(' ')[1])
                errXCentre = float(lines[2].split(':')[1].split(' ')[1])
                errYCentre = float(lines[4].split(':')[1].split(' ')[1])
        continue
        # get polar histogram:
        if not os.path.exists(eventPath + "/histoPolarWindow.png"):
            (polarHistoWindow) = getPolarHisto(eventPath, configFile, "histoPolarWindow", photonCounterDataWindow, bestXCircleCenter, bestYCircleCenter, startWindow, endWindow)
        #if not os.path.exists(eventPath + "/histoPolarELVE.png"):
            #(polarHistoFrame) = getPolarHisto(eventPath, configFile, "histoPolarELVE", photonCounterDataWindow, bestXCircleCenter, bestYCircleCenter, startWindow, configFile.FrameEnd)

        # update pixelsInfos:
        if not os.path.exists(eventPath + "/infoPixelPhi.txt"):
            getPhiEnergy(eventPath, pixelsInfos, bestXCircleCenter, bestYCircleCenter)

        #search for second ring:
        #pixelsInfosII = getPixelsInfosII(eventPath + "/pixelsII", configFile, photonCounterDataWindow, pixelsInfos, startWindow)
