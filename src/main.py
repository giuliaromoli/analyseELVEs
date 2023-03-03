# created:15/02/2023

import ROOT
import os

from src.utils.getWindowData import getEventInfos
from src.utils.getPixelsInfos import getPixelsInfos
from src.utils.getPixelsInfos import readPixelsInfos
from src.utils.getPixelsInfos import getPixelsPlots
from src.utils.getCircleCenter import getCircleCenter
from src.utils.fillHisto import getPolarHisto

from easydict import EasyDict as ed
import yaml
from src.utils.plotELVEResults import getPhiEnergy
from src.utils.plotELVEResults import readPixelsInfosPhi
from src.utils.plotELVEResults import plotResults
from src.utils.plotELVEResults import subtractBackground

ROOT.gROOT.SetBatch(ROOT.kTRUE)  # do not open figures!

if __name__ == '__main__':

    eventLabel = "30-12-2019"
    print('start analysing ' + eventLabel)

    eventPath = r'/home/jule/Scrivania/Mini-Euso/analyseELVEs/results/' + eventLabel + '/'
    if not os.path.exists(eventPath):
        os.makedirs(eventPath)

    ### GET PIXELS DATA:

    if not os.path.exists(eventPath + "pixels/infoPixels.txt"):
        (photonCounterDataWindow, configFile, startWindow, endWindow) = getEventInfos(eventLabel, eventPath)
        # fit pixels with GaussExp:
        pixelsInfos = getPixelsInfos(eventPath + "pixels", configFile, photonCounterDataWindow, startWindow)
        getPixelsPlots(eventPath + "/pixelsPlots", pixelsInfos, configFile)
        # fit circles:
        (bestXCircleCenter, bestYCircleCenter, errCentreX, errCentreY) = getCircleCenter(
            eventPath + "circleFrames/", pixelsInfos, photonCounterDataWindow, configFile, startWindow)
        # get polar histograms:
        if not os.path.exists(eventPath + "/histoPolarWindow.png"):
            (polarHistoWindow) = getPolarHisto(eventPath, configFile, "histoPolarWindow", photonCounterDataWindow,
                                               bestXCircleCenter, bestYCircleCenter, startWindow, endWindow)
        if not os.path.exists(eventPath + "/histoPolarELVE.png"):
            (polarHistoFrame) = getPolarHisto(eventPath, configFile, "histoPolarELVE", photonCounterDataWindow,
                                              bestXCircleCenter, bestYCircleCenter, startWindow, configFile.FrameEnd)
    else:
        pixelsInfos = readPixelsInfos(eventPath)

    ### PLOT ELVE RESULTS:

    # get config file:
    with open(
            r'/home/jule/Scrivania/Mini-Euso/analyseELVEs/src/configResults/configResults-' + eventLabel + '.yaml') as file:
        configR = yaml.safe_load(file)
    configRFile = ed(configR)

    # get infos:
    if not os.path.exists(eventPath + "infoPixelPhi.txt"):
        # get distance and angle
        getPhiEnergy(eventPath, pixelsInfos, configRFile.CircleCoords.XCoord, configRFile.CircleCoords.YCoord)
        pixelsInfosPhi = readPixelsInfosPhi(eventPath)
    else:
        pixelsInfosPhi = readPixelsInfosPhi(eventPath)

    # plot results:
    plotResults(pixelsInfosPhi, configRFile, eventPath + "resultsWithBackground/")
    # plot results (with background subtraction):
    pixelsInfosPhiSub = subtractBackground(pixelsInfosPhi)
    plotResults(pixelsInfosPhiSub, configRFile, eventPath + "resultsWithoutBackground/")
