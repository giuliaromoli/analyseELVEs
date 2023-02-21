import ROOT
import numpy as np
import csv
import math


def getPixelsPositions(pixelFileName, HH):
    kmdistanceX = np.zeros((48, 48))  # pixel centre
    kmdistanceY = np.zeros((48, 48))
    with open(r'/home/jule/Scrivania/Mini-Euso/analyseELVEs/src/utils/' + pixelFileName + '.txt') as PixelFile:
        reader = csv.reader(PixelFile)
        for row in reader:  # for each pixel
            pixX = int(row[0])  # pixel X coordinate
            pixY = int(row[1])  # pixel Y coordinate
            th0 = float(row[2])  # Theta angle pixel centre
            ph0 = float(row[3])  # Phi angle pixel center
            kmdistanceX[pixX][pixY] = HH * math.tan(th0) * math.cos(ph0)
            kmdistanceY[pixX][pixY] = HH * math.tan(th0) * math.sin(ph0)
    return kmdistanceX, kmdistanceY


def fillHisto(newPath, pixelFileName, HH, infoHisto, minHisto, maxHisto, markerPosX, markerPosY, color, title, figureName):
    kmdistanceX = np.zeros((48, 48))  # pixel centre
    kmdistanceY = np.zeros((48, 48))
    kmdistanceX1 = np.zeros((48, 48))  # RightTop (ETOS view)
    kmdistanceY1 = np.zeros((48, 48))
    kmdistanceX2 = np.zeros((48, 48))  # LeftTop
    kmdistanceY2 = np.zeros((48, 48))
    kmdistanceX3 = np.zeros((48, 48))  # LeftBottom
    kmdistanceY3 = np.zeros((48, 48))
    kmdistanceX4 = np.zeros((48, 48))  # RightBottom
    kmdistanceY4 = np.zeros((48, 48))

    with open(r'/home/jule/Scrivania/Mini-Euso/analyseELVEs/src/utils/' + pixelFileName + '.txt') as PixelFile:
        reader = csv.reader(PixelFile)
        for row in reader:  # for each pixel

            pixX = int(row[0])  # pixel X coordinate
            pixY = int(row[1])  # pixel Y coordinate
            th0 = float(row[2])  # Theta angle pixel centre
            ph0 = float(row[3])  # Phi angle pixel center
            th1 = float(row[4])  # Theta angle right-top vertex
            ph1 = float(row[5])  # Phi angle right-top vertex
            th2 = float(row[6])
            ph2 = float(row[7])
            th3 = float(row[8])
            ph3 = float(row[9])
            th4 = float(row[10])
            ph4 = float(row[11])

            kmdistanceX[pixX][pixY] = HH * math.tan(th0) * math.cos(ph0)
            kmdistanceY[pixX][pixY] = HH * math.tan(th0) * math.sin(ph0)
            kmdistanceX1[pixX][pixY] = HH * math.tan(th1) * math.cos(ph1)
            kmdistanceY1[pixX][pixY] = HH * math.tan(th1) * math.sin(ph1)
            kmdistanceX2[pixX][pixY] = HH * math.tan(th2) * math.cos(ph2)
            kmdistanceY2[pixX][pixY] = HH * math.tan(th2) * math.sin(ph2)
            kmdistanceX3[pixX][pixY] = HH * math.tan(th3) * math.cos(ph3)
            kmdistanceY3[pixX][pixY] = HH * math.tan(th3) * math.sin(ph3)
            kmdistanceX4[pixX][pixY] = HH * math.tan(th4) * math.cos(ph4)
            kmdistanceY4[pixX][pixY] = HH * math.tan(th4) * math.sin(ph4)

    # plot timeDelay/pixel:
    canvas = ROOT.TCanvas()
    histo = ROOT.TH2Poly()
    histo.SetMaximum(maxHisto)
    histo.SetMinimum(minHisto)
    for px in range(48):
        for py in range(48):
            tempXvec = np.zeros(4)
            tempYvec = np.zeros(4)
            tempXvec[0] = kmdistanceX1[px][py]
            tempXvec[1] = kmdistanceX2[px][py]
            tempXvec[2] = kmdistanceX3[px][py]
            tempXvec[3] = kmdistanceX4[px][py]
            tempYvec[0] = kmdistanceY1[px][py]
            tempYvec[1] = kmdistanceY2[px][py]
            tempYvec[2] = kmdistanceY3[px][py]
            tempYvec[3] = kmdistanceY4[px][py]
            histo.AddBin(4, tempXvec, tempYvec)
    for px in range(48):
        for py in range(48):
            histo.Fill(-kmdistanceX[px][py], -kmdistanceY[px][py], infoHisto[px][py])
    histo.SetTitle(title)
    histo.Draw("COLZ")
    histo.SetStats(0)
    for k in range(len(markerPosX)):
        m = ROOT.TMarker()
        m.SetMarkerStyle(20)
        m.SetMarkerColor(color)
        m.DrawMarker(markerPosX[k], markerPosY[k])
    canvas.Update()
    canvas.SetCanvasSize(800, 800)
    canvas.SaveAs(newPath + figureName)
    return histo, canvas


def getPolarHisto(newPath, configFile, namePlot, data, bestCentreX, bestCentreY, start, end):

    binEnd = configFile.NumBinsToPlot * configFile.BinSize + configFile.BinStart
    plotHistoPolar = ROOT.TCanvas("plotHistoPolar")
    polarLim = end - start + 1
    histoPolar = ROOT.TH2F("dueddipolar", "Polar Histogram; time [#frame]; radius [km]", int(polarLim), start, end,
                           configFile.BumBinsToPlot, configFile.BinStart, binEnd)

    (kmdistanceX, kmdistanceY) = getPixelsPositions(configFile.PixelFileName, configFile.HH)

    for time in range(start, end + 1):
        for iix in range(48):
            for iiy in range(48):
                distanceElfePixel = np.sqrt(
                    np.power(-kmdistanceX[iix][iiy] - bestCentreX, 2) + np.power(-kmdistanceY[iix][iiy] - bestCentreY,
                                                                                 2))
                histoPolar.Fill(time, distanceElfePixel, data[0][0][time - start][iix][iiy])
    histoPolar.Draw("COLZ")
    histoPolar.SetStats(0)
    # histoPolar.GetZaxis().SetRangeUser(0,100)
    plotHistoPolar.SaveAs(newPath + "/" + namePlot + ".png")
    return histoPolar


if __name__ == '__main__':
    pass
