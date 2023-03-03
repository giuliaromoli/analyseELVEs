import ROOT
import csv
import math
from easydict import EasyDict as ed
import numpy as np
import os
from src.utils.fillHisto import fillHisto
from src.utils.fillHisto import getPixelsPositions
import yaml


def getEventInfos(eventLabel, eventPath):

    # get config file:
    with open(r'/home/jule/Scrivania/Mini-Euso/analyseELVEs/src/config/config-' + eventLabel + '.yaml') as file:
        config = yaml.safe_load(file)
    configFile = ed(config)

    # get window limits:
    if not os.path.exists(eventPath + "infoWindow/infoWindow.txt"):
        (startWindow, endWindow) = getWindowLimits(eventPath, configFile)
    else:
        with open(eventPath + "infoWindow/infoWindow.txt") as f:
            lines = f.readlines()
            startWindow = int(lines[0].split(':')[1])
            endWindow = int(lines[1].split(':')[1])
    # print(startWindow,endWindow)

    # get data:
    photonCounterDataWindow = getData(startWindow, endWindow, eventPath, configFile)

    return photonCounterDataWindow, configFile, startWindow, endWindow


def getPixelsDeltaTime(eventPath, pixelFileName, applyTimeCorr, HH):
    deltaTime = np.zeros((48, 48))  # time delay
    deltaTimeFrame = np.zeros((48, 48))
    with open(r'/home/jule/Scrivania/Mini-Euso/analyseELVEs/src/utils/' + pixelFileName + '.txt') as PixelFile:
        reader = csv.reader(PixelFile)
        for row in reader:  # for each pixel
            pixX = int(row[0])  # pixel X coordinate
            pixY = int(row[1])  # pixel Y coordinate
            th0 = float(row[2])  # Theta angle pixel centre
            # get time delay with respect to center of focal surface:
            if applyTimeCorr == 1:
                deltaTime[pixX][pixY] = HH / 3e5 * (1 / np.cos(th0) - 1)
                deltaTimeFrame[pixX][pixY] = deltaTime[pixX][pixY] * 1e6 / 2.5
            else:
                deltaTime[pixX][pixY] = 0
                deltaTimeFrame[pixX][pixY] = 0
    # plot timeDelay/pixel:
    if not os.path.exists(eventPath + "/deltaTime.png"):
        (histo, canvas) = fillHisto(eventPath, pixelFileName, HH, deltaTimeFrame, 0, int(np.max(deltaTimeFrame)), [],
                                    [],
                                    0, "delta time [GTU]", "/deltaTime.png")
    return deltaTimeFrame


def getWindowLimits(eventPath, configFile):
    # load .root file name:
    rootFileName = configFile.RootFileName
    version = configFile.Version
    fullRootFileName = rootFileName + version
    inputFile = r'/home/jule/Scrivania/Mini-Euso/analyseELVEs/inputFiles/' + fullRootFileName + '.root'
    # print(inputFile)
    # open .root input file:
    rootFile = ROOT.TFile(inputFile, "read")
    # aquire data (D1 sampling):
    tevent = rootFile.Get("tevent");
    nentries = tevent.GetEntries()  # number of entries
    photonCounterData = np.zeros((1, 1, 48, 48), dtype=np.float32)  # photon data (uncorr)
    tevent.SetBranchAddress("photon_count_data", photonCounterData)  # set address
    # open .txt file:
    newPath = eventPath + "/infoWindow"
    if not os.path.exists(newPath):
        os.makedirs(newPath)
    fileInfoWindow = open(newPath + "/infoWindow.txt", "w+")
    # define temporal limits [GTUs]:
    frameStart = configFile.FrameStart
    arrayDiffUnixTime = np.ndarray(nentries)
    arrayUnixTime = np.ndarray(nentries)
    diffUnixTimeCanvas = ROOT.TCanvas("diffUnix")
    diffUnixTimePlot = ROOT.TGraph()
    for t in range(nentries):
        tevent.GetEntry(t)
        arrayUnixTime[t] = tevent.timestamp_unix
        if t > 0:
            arrayDiffUnixTime[t] = arrayUnixTime[t] - arrayUnixTime[t - 1]
            diffUnixTimePlot.SetPoint(diffUnixTimePlot.GetN(), t, arrayDiffUnixTime[t])
    diffUnixTimePlot.SetMarkerStyle(20)
    diffUnixTimePlot.SetMarkerColor(4)
    diffUnixTimePlot.SetMarkerSize(1)
    diffUnixTimePlot.Draw("AP")
    windowsStarts = np.asarray(np.where(arrayDiffUnixTime > 0))
    startFullWindow = windowsStarts[0][max(max(np.where(windowsStarts[0][:] - frameStart < 0)))]
    endFullWindow = windowsStarts[0][max(max(np.where(windowsStarts[0][:] - frameStart < 0))) + 1]
    prevWindow = windowsStarts[0][max(max(np.where(windowsStarts[0][:] - frameStart < 0))) - 1]
    unixStart = arrayUnixTime[startFullWindow]
    unixEnd = arrayUnixTime[endFullWindow]
    unixPrev = arrayUnixTime[prevWindow]
    nTriggers = (endFullWindow - startFullWindow) / 128
    fileInfoWindow.write("Start window [GTU]: %d \nEnd window [GTU]: %d \n" % (startFullWindow, endFullWindow))
    fileInfoWindow.write("GTUs in window: %d \n" % (endFullWindow - startFullWindow))
    fileInfoWindow.write("# triggers: %d \n" % nTriggers)
    fileInfoWindow.write("window unix time [s]: %d \n" % unixStart)
    fileInfoWindow.write("seconds from prev window: %d \n" % (unixStart - unixPrev))
    fileInfoWindow.write("seconds to next window: %d \n" % (unixEnd - unixStart))
    lineStart = ROOT.TLine(startFullWindow, 0, startFullWindow, 1e10);
    lineStart.SetLineColor(4)
    lineStart.SetLineStyle(9)
    lineStart.Draw()
    lineEnd = ROOT.TLine(endFullWindow, 0, endFullWindow, 1e10);
    lineEnd.SetLineColor(4)
    lineEnd.SetLineStyle(9)
    lineEnd.Draw()
    diffUnixTimePlot.GetXaxis().SetTitle("time [D1 GTU]")
    diffUnixTimePlot.GetYaxis().SetNoExponent(1)
    diffUnixTimePlot.GetYaxis().SetTitle("diff unix time [s]")
    diffUnixTimeCanvas.SaveAs(newPath + "/" + "diffUnix.png")
    # plot array unix times:
    unixTimeCanvas = ROOT.TCanvas("unix")
    unixTimePlot = ROOT.TGraph()
    for t in range(nentries):
        unixTimePlot.SetPoint(unixTimePlot.GetN(), t, arrayUnixTime[t])
    unixTimePlot.SetMarkerStyle(20)
    unixTimePlot.SetMarkerColor(4)
    unixTimePlot.SetMarkerSize(1)
    unixTimePlot.Draw()
    lineStart = ROOT.TLine(startFullWindow, 0, startFullWindow, 1e10);
    lineStart.SetLineColor(4)
    lineStart.SetLineStyle(9)
    lineStart.Draw()
    lineEnd = ROOT.TLine(endFullWindow, 0, endFullWindow, 1e10);
    lineEnd.SetLineColor(4)
    lineEnd.SetLineStyle(9)
    lineEnd.Draw()
    unixTimePlot.GetXaxis().SetTitle("time [D1 GTU]")
    unixTimePlot.GetYaxis().SetNoExponent(1)
    unixTimePlot.GetYaxis().SetTitle("unix time [s]")
    unixTimeCanvas.SaveAs(newPath + "/" + "unix.png")
    # get trigger frames:
    for k in range(4):
        for time in ([startFullWindow + (128 * (k + 1) - 1), startFullWindow + (128 * (k + 1)),
                      startFullWindow + (128 * (k + 1)) + 1]):
            tevent.GetEntry(time)
            (histo, canvas) = fillHisto(newPath, configFile.PixelFileName, configFile.HH, photonCounterData[0][0][:][:],
                                        0, configFile.MaxHisto, [], [], 0,
                                        "frame " + str(time), "/frame" + str(time) + ".png")

    fileInfoWindow.close()
    return (startFullWindow, endFullWindow)


def getData(startWindow, endWindow, eventPath, configFile):
    # load .root file name:
    rootFileName = configFile.RootFileName
    version = configFile.Version
    fullRootFileName = rootFileName + version
    inputFile = r'/home/jule/Scrivania/Mini-Euso/analyseELVEs/inputFiles/' + fullRootFileName + '.root'
    # print(inputFile)
    # open .root input file:
    rootFile = ROOT.TFile(inputFile, "read")
    # acquire data (D1 sampling):
    tevent = rootFile.Get("tevent")
    photonCounterData = np.zeros((1, 1, 48, 48), dtype=np.float32)
    tevent.SetBranchAddress("photon_count_data", photonCounterData)
    # apply time correction:
    deltaTimeFrame = getPixelsDeltaTime(eventPath, configFile.PixelFileName, configFile.ApplyTimeCorr, configFile.HH)
    DELTAMAX = int(np.max(deltaTimeFrame))  # max time delay
    photonCounterData3D = np.zeros((1, 1, endWindow + 1 - startWindow + DELTAMAX, 48, 48), dtype=np.float32)
    n = tevent.GetEntries()
    tempCounter = DELTAMAX
    for x in range(n):
        if startWindow <= x <= endWindow:
            tevent.GetEntry(x)
            for iix in range(48):
                for iiy in range(48):
                    photonCounterData3D[0][0][tempCounter - int(deltaTimeFrame[iix][iiy])][iix][iiy] = \
                        photonCounterData[0][0][iix][iiy]
            tempCounter += 1
    # get data:
    photonCounterDataWindow = np.zeros((1, 1, endWindow + 1 - startWindow, 48, 48), dtype=np.float32)
    for time in range(startWindow, endWindow + 1):
        photonCounterDataWindow[0, 0, time - startWindow, :, :] = photonCounterData3D[0, 0, time - startWindow, :, :]
        photonCounterDataWindow[photonCounterDataWindow < 0] = 0
    # apply moving average:
    jump = 0
    totLen = len(photonCounterDataWindow[0, 0, :, :, :])
    photonCounterDataWindowAvg = np.zeros((1, 1, int(totLen), 48, 48), dtype=np.float32)
    for px in range(48):
        for py in range(48):
            for time in range(jump, totLen - (jump + 1)):
                photonCounterDataWindowAvg[0, 0, time, px, py] = np.average(
                    photonCounterDataWindow[0, 0, time - jump:time + jump + 1, px,
                    py])  # ,weights=np.array([0.155,0.69,0.155])

    # plot pixels' signals:
    newPath = eventPath + "/pixelsProfiles"
    if not os.path.exists(newPath):
        os.makedirs(newPath)
    (kmdistanceX, kmdistanceY) = getPixelsPositions(configFile.PixelFileName, configFile.HH)
    for px in range(48):
        for py in range(48):
            pixelLabel = px + 48 * py
            if not os.path.exists(newPath + "/" + str(pixelLabel) + ".png"):
                xCoord = -kmdistanceX[px][py]  # x-coord
                yCoord = -kmdistanceY[px][py]  # y-coord
                maxEnergy = max(photonCounterDataWindowAvg[0, 0, :, px, py])
                if ((configFile.limXmin <= xCoord <= configFile.limXmax) and (
                        configFile.limYmin <= yCoord <= configFile.limYmax) and (maxEnergy > 0)):
                    plotPIXProfile = ROOT.TCanvas("plotPixProfile" + str(pixelLabel))
                    plotPIX = ROOT.TMultiGraph()
                    plotEnePixel = ROOT.TGraph()
                    for time in range(len(photonCounterDataWindowAvg[0, 0, :, px, py])):
                        plotEnePixel.SetPoint(plotEnePixel.GetN(), time + startWindow,
                                              photonCounterDataWindowAvg[0, 0, time, px, py])
                    plotEnePixel.SetMarkerStyle(20)
                    plotEnePixel.SetMarkerColor(1)
                    plotEnePixel.SetMarkerSize(0.5)
                    plotPIX.Add(plotEnePixel)
                    plotPIX.Draw("APL")
                    plotPIXProfile.SaveAs(newPath + "/" + str(pixelLabel) + ".png")
    return photonCounterDataWindow


if __name__ == '__main__':
    pass
