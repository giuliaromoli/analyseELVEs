import ROOT
import numpy as np
import math
import os
from array import array
import scipy.optimize as sciopt
from src.utils.fillHisto import fillHisto
from scipy.optimize import curve_fit
import scipy.special as sse


def getResults(pixelsInfosPhi, configRFile, eventPath):
    # plot radius results:
    plotRadiusResults(pixelsInfosPhi, configRFile, eventPath + "radius/")
    # plot energy results (without background subtraction):
    plotResults(pixelsInfosPhi, configRFile, eventPath + "energy/withoutBackground/")
    # plot energy results (with background subtraction):
    pixelsInfosPhiSub = subtractBackground(pixelsInfosPhi)
    plotResults(pixelsInfosPhiSub, configRFile, eventPath + "energy/withBackground/")


def subtractBackground(pixelsInfosPhi):
    pixelsInfosPhiSub = []
    for ll in range(len(pixelsInfosPhi)):
        pixelsInfosPhiSub.append(list(pixelsInfosPhi[ll]))
    for ll in range(len(pixelsInfosPhiSub)):
        pixelsInfosPhiSub[ll][7] = pixelsInfosPhiSub[ll][7] - pixelsInfosPhiSub[ll][25]
        pixelsInfosPhiSub[ll][11] = pixelsInfosPhiSub[ll][11] - pixelsInfosPhiSub[ll][25]
        pixelsInfosPhiSub[ll][13] = pixelsInfosPhiSub[ll][13] - pixelsInfosPhiSub[ll][25]
    return pixelsInfosPhiSub


def getPhiEnergy(newPath, pixelsInfos, ElveCentreX, ElveCentreY):
    pixelsInfosPhi = ()
    fileInfoPixelPhi = open(newPath + "/infoPixelPhi.txt", "w+")  # open debugFile

    for t in range(len(pixelsInfos)):
        xCoord = pixelsInfos[t][3] - ElveCentreX
        yCoord = pixelsInfos[t][4] - ElveCentreY
        distanceFromCenter = np.sqrt(np.power(xCoord, 2) + np.power(yCoord, 2))
        thetaRad = np.arctan(yCoord / xCoord)
        thetaDeg = thetaRad * 180 / math.pi
        if xCoord < 0:
            thetaDeg = thetaDeg + 180
        pixelsInfosPhi.append((pixelsInfos[t][0], pixelsInfos[t][1], pixelsInfos[t][2], pixelsInfos[t][3],
                               pixelsInfos[t][4], pixelsInfos[t][5], pixelsInfos[t][6], pixelsInfos[t][7],
                               pixelsInfos[t][8], pixelsInfos[t][9], pixelsInfos[t][10], pixelsInfos[t][11],
                               pixelsInfos[t][12], pixelsInfos[t][13], pixelsInfos[t][14], pixelsInfos[t][15],
                               pixelsInfos[t][16], pixelsInfos[t][17], pixelsInfos[t][18], pixelsInfos[t][19],
                               pixelsInfos[t][20], pixelsInfos[t][21], pixelsInfos[t][22], distanceFromCenter, thetaDeg,
                               pixelsInfos[t][23], pixelsInfos[t][24], pixelsInfos[t][25], pixelsInfos[t][26]))
        fileInfoPixelPhi.write(
            "%d \t %d \t %d \t %5.2f \t %5.2f \t %d \t %d \t %5.2f \t %5.2f \t %d \t %d \t %5.2f \t %5.2f \t %5.2f \t %5.2f \t %5.2f \t %5.2f \t %d \t %5.2f \t %5.2f \t %5.2f \t %5.2f \t %5.2f \t %5.2f \t %5.2f \t %5.2f \t %5.2f \t %d \t %d \n" % (
                pixelsInfos[t][0], pixelsInfos[t][1], pixelsInfos[t][2], pixelsInfos[t][3], pixelsInfos[t][4],
                pixelsInfos[t][5], pixelsInfos[t][6], pixelsInfos[t][7], pixelsInfos[t][8], pixelsInfos[t][9],
                pixelsInfos[t][10], pixelsInfos[t][11], pixelsInfos[t][12], pixelsInfos[t][13], pixelsInfos[t][14],
                pixelsInfos[t][15], pixelsInfos[t][16], pixelsInfos[t][17], pixelsInfos[t][18], pixelsInfos[t][19],
                pixelsInfos[t][20], pixelsInfos[t][21], pixelsInfos[t][22], distanceFromCenter, thetaDeg,
                pixelsInfos[t][23], pixelsInfos[t][24], pixelsInfos[t][25], pixelsInfos[t][26]))
    fileInfoPixelPhi.close()
    return pixelsInfosPhi


def readPixelsInfosPhi(eventPath):
    pixelsInfosPhi = []
    if not os.path.exists(eventPath + "infoPixelPhi.txt"):
        print("no pixelsInfoPhi.txt file found!")
        return pixelsInfosPhi
    else:
        with open(eventPath + "/infoPixelPhi.txt") as f:
            for line in f:
                pixelLabel = int(line.split('\t')[0])
                pixX = int(line.split('\t')[1])
                pixY = int(line.split('\t')[2])
                x = float(line.split('\t')[3])
                y = float(line.split('\t')[4])
                nPeak = int(line.split('\t')[5])
                timeMax = int(line.split('\t')[6])
                maxEnergy = float(line.split('\t')[7])
                sigmaMaxEnergy = float(line.split('\t')[8])
                fullWidth = int(line.split('\t')[9])
                deExcTime = int(line.split('\t')[10])
                totEnergy = float(line.split('\t')[11])
                sigmaTotEnergy = float(line.split('\t')[12])
                energyPeak = float(line.split('\t')[13])
                sigmaEnergyPeak = float(line.split('\t')[14])
                # from fit with gaussExp:
                tauFit = float(line.split('\t')[15])
                sigmaFit = float(line.split('\t')[16])
                meanFit = int(line.split('\t')[17])
                areaFit = float(line.split('\t')[18])
                maxEnergyFit = float(line.split('\t')[19])
                sigmaMaxEnergyFit = float(line.split('\t')[20])
                energyPeakFit = float(line.split('\t')[21])
                sigmaEnergyPeakFit = float(line.split('\t')[22])
                distanceFromCenter = float(line.split('\t')[23])
                thetaAngle = float(line.split('\t')[24])
                thresholdSinglePeak = float(line.split('\t')[25])
                lambdaFit = float(line.split('\t')[26])
                totDuration = int(line.split('\t')[27])
                totDeExcTime = int(line.split('\t')[28].split('\n')[0])
                pixelsInfosPhi.append((pixelLabel, pixX, pixY, x, y, nPeak, timeMax, maxEnergy, sigmaMaxEnergy,
                                       fullWidth, deExcTime, totEnergy, sigmaTotEnergy, energyPeak, sigmaEnergyPeak,
                                       tauFit, sigmaFit, meanFit, areaFit, maxEnergyFit, sigmaMaxEnergyFit,
                                       energyPeakFit, sigmaEnergyPeakFit, distanceFromCenter, thetaAngle,
                                       thresholdSinglePeak, lambdaFit, totDuration, totDeExcTime))
    return pixelsInfosPhi


def plotPixelsResults(pixelsInfosPhi, configRFile, newPath):
    if not os.path.exists(newPath):
        os.makedirs(newPath)
    fileInfoFit = open(newPath + "/infoFits.txt", "w+")

    ### SELECT RING:
    selPixelsInfosPhi = []
    for ll in range(len(pixelsInfosPhi)):
        if pixelsInfosPhi[ll][5] == 1:
            selPixelsInfosPhi.append(pixelsInfosPhi[ll])

    ### DURATION:
    canvas = ROOT.TCanvas("duration")
    multiGraph = plotResult(selPixelsInfosPhi, 6, "time [GTU]", configRFile.Frames.Start, configRFile.Frames.End, 1, 27,
                            "duration [GTU]", 0, 15, "both", "none", "none", configRFile, fileInfoFit)
    canvas.SaveAs(newPath + "duration.png")

    ### DE-EXCITATION:


def saveFigure(canvas, path, name):
    canvas.SaveAs(path + name + ".png")
    canvas.SaveAs(path + name + ".root")


def plotRadiusResults(pixelsInfosPhi, configRFile, newPath):

    if not os.path.exists(newPath):
        os.makedirs(newPath)
    fileInfoFit = open(newPath + "/infoFits.txt", "w+")

    ### SELECT RING:
    selPixelsInfosPhi = []
    for ll in range(len(pixelsInfosPhi)):
        if pixelsInfosPhi[ll][5] == 1:
            selPixelsInfosPhi.append(pixelsInfosPhi[ll])

    ### RADIUS:
    fileInfoFit.write("---RADIUS VS TIME \n")
    canvas = ROOT.TCanvas("radius")
    legend = ROOT.TLegend(0.1, 0.7, 0.45, 0.9)
    multiGraph = plotResult(selPixelsInfosPhi, 6, "time [GTU]", configRFile.Radius.XMin, configRFile.Radius.XMax, configRFile.Radius.Jump, 23,
                            "distance [km]", configRFile.Radius.YMin, configRFile.Radius.YMax, 0, 360, "both", "radius",
                            "scatter", configRFile, fileInfoFit, legend, 20, 1)
    legend.SetTextSize(0.04)
    legend.Draw()
    saveFigure(canvas, newPath, "radius")
    fileInfoFit.write("\n")
    fileInfoFit.close()

    ### PIXELS VS TIME:
    canvas = ROOT.TCanvas("pixels")
    minTime = configRFile.Radius.XMin
    maxTime = configRFile.Radius.XMax
    histoPixels = ROOT.TH1F("TD", "", maxTime - minTime + 1, minTime, maxTime)
    for t in range(len(selPixelsInfosPhi)):
        histoPixels.Fill(selPixelsInfosPhi[t][6], 1)
    histoPixels.SetStats(0)
    histoPixels.SetFillColor(38)
    histoPixels.SetFillStyle(3001)
    histoPixels.Draw()
    histoPixels.GetXaxis().SetTitle("time [GTU]")
    histoPixels.GetYaxis().SetTitle("# pixels")
    saveFigure(canvas, newPath, "pixelsVsTime")

    ### PIXELS VS DISTANCE:
    canvas = ROOT.TCanvas("pixels")
    minTime = configRFile.Radius.YMin
    maxTime = configRFile.Radius.YMax
    histoPixels = ROOT.TH1F("TD", "", maxTime - minTime + 1, minTime, maxTime)
    for t in range(len(selPixelsInfosPhi)):
        histoPixels.Fill(selPixelsInfosPhi[t][23], 1)
    histoPixels.SetStats(0)
    histoPixels.SetFillColor(38)
    histoPixels.SetFillStyle(3001)
    histoPixels.Draw()
    histoPixels.GetXaxis().SetTitle("distance [km]")
    histoPixels.GetYaxis().SetTitle("# pixels")
    saveFigure(canvas, newPath, "pixelsVsDistance")


def plotEnergyResultsVsTime(selPixelsInfosPhi, configRFile, newPath):

    if not os.path.exists(newPath):
        os.makedirs(newPath)
    fileInfoFit = open(newPath + "/infoFits.txt", "w+")

    fileInfoFit.write("---MAX ENERGY VS TIME \n")
    canvas = ROOT.TCanvas("energyMax")
    legendT = ROOT.TLegend(0.65, 0.55, 0.9, 0.9)
    multiGraph = plotResult(selPixelsInfosPhi, 6, "time [GTU]", configRFile.EnergyTime.XMin,
                            configRFile.EnergyTime.XMax,
                            configRFile.EnergyTime.Jump, 7,
                            "energy max [ADC]", configRFile.EnergyTime.YMin, configRFile.EnergyTime.YMax, 0, 360,
                            "both",
                            "energyTime", "bar", configRFile, fileInfoFit, legendT, 20, 1)
    legendT.SetTextSize(0.04)
    legendT.Draw()
    saveFigure(canvas, newPath, "energyMaxVsTime" + str(configRFile.EnergyTime.Jump))
    fileInfoFit.write("\n")

    # in log scale:
    canvas = ROOT.TCanvas("energyMaxLog")
    gPad = ROOT.TPad("pad1", "pad1", 0.1, 0.1, 0.9, 0.9, 0)
    gPad.Draw()
    gPad.cd()
    gPad.SetLogy()
    legendTL = ROOT.TLegend(0.7, 0.8, 0.9, 0.9)
    multiGraph = plotResult(selPixelsInfosPhi, 6, "time [GTU]", configRFile.EnergyTime.XMin,
                            configRFile.EnergyTime.XMax,
                            configRFile.EnergyTime.Jump, 7,
                            "energy max [ADC]", configRFile.EnergyTime.YMin, configRFile.EnergyTime.YMax, 0, 360,
                            "both",
                            "energyTimeLog", "bar", configRFile, fileInfoFit, legendTL, 20, 1)
    legendTL.SetTextSize(0.045)
    legendTL.Draw()
    saveFigure(canvas, newPath, "energyMaxVsTime" + str(configRFile.EnergyTime.Jump) + "Log")
    fileInfoFit.write("\n")

    fileInfoFit.write("---PEAK ENERGY VS TIME \n")
    canvas = ROOT.TCanvas("energyPeak")
    legendP = ROOT.TLegend(0.65, 0.55, 0.9, 0.9)
    multiGraph = plotResult(selPixelsInfosPhi, 6, "time [GTU]", configRFile.EnergyTime.XMin,
                            configRFile.EnergyTime.XMax,
                            configRFile.EnergyTime.Jump, 13,
                            "energy peak [ADC]", configRFile.EnergyTime.YMin, configRFile.EnergyTime.YMax, 0, 360,
                            "both",
                            "energyTime", "bar", configRFile, fileInfoFit, legendP, 20, 1)
    legendP.SetTextSize(0.04)
    legendP.Draw()
    saveFigure(canvas, newPath, "energyPeakVsTime" + str(configRFile.EnergyTime.Jump))
    fileInfoFit.write("\n")

    # in log scale:
    canvas = ROOT.TCanvas("energyPeakLog")
    gPad = ROOT.TPad("pad1", "pad1", 0.1, 0.1, 0.9, 0.9, 0)
    gPad.Draw()
    gPad.cd()
    gPad.SetLogy()
    legendPL = ROOT.TLegend(0.7, 0.8, 0.9, 0.9)
    multiGraph = plotResult(selPixelsInfosPhi, 6, "time [GTU]", configRFile.EnergyTime.XMin,
                            configRFile.EnergyTime.XMax,
                            configRFile.EnergyTime.Jump, 13,
                            "energy peak [ADC]", configRFile.EnergyTime.YMin, configRFile.EnergyTime.YMax, 0, 360,
                            "both",
                            "energyTimeLog", "bar", configRFile, fileInfoFit, legendPL, 20, 1)
    legendPL.SetTextSize(0.045)
    legendPL.Draw()
    saveFigure(canvas, newPath, "energyPeakVsTime" + str(configRFile.EnergyTime.Jump) + "Log")
    fileInfoFit.write("\n")
    fileInfoFit.close()


def plotEnergyResultsVsDistance(selPixelsInfosPhi, configRFile, newPath):

    if not os.path.exists(newPath):
        os.makedirs(newPath)
    fileInfoFit = open(newPath + "/infoFits.txt", "w+")

    fileInfoFit.write("---MAX ENERGY VS DISTANCE \n")
    canvas = ROOT.TCanvas("energyRMax")
    legendR = ROOT.TLegend(0.75, 0.55, 0.9, 0.9)
    multiGraph = plotResult(selPixelsInfosPhi, 23, "distance [km]", configRFile.EnergyRadius.XMin,
                            configRFile.EnergyRadius.XMax,
                            configRFile.EnergyRadius.Jump, 7,
                            "energy max [ADC]", configRFile.EnergyRadius.YMin, configRFile.EnergyRadius.YMax, 0, 360,
                            "both",
                            "energyRadius", "bar", configRFile, fileInfoFit, legendR, 20, 1)
    legendR.SetTextSize(0.04)
    legendR.Draw()
    saveFigure(canvas, newPath, "energyMaxVsDistance" + str(configRFile.EnergyTime.Jump))
    fileInfoFit.write("\n")

    # in log scale:
    canvas = ROOT.TCanvas("energyRMaxLog")
    gPad = ROOT.TPad("pad1", "pad1", 0.1, 0.1, 0.9, 0.9, 0)
    gPad.Draw()
    gPad.cd()
    gPad.SetLogy()
    legendRL = ROOT.TLegend(0.7, 0.8, 0.9, 0.9)
    multiGraph = plotResult(selPixelsInfosPhi, 23, "distance [km]", configRFile.EnergyRadius.XMin,
                            configRFile.EnergyRadius.XMax,
                            configRFile.EnergyRadius.Jump, 7,
                            "energy max [ADC]", configRFile.EnergyRadius.YMin, configRFile.EnergyRadius.YMax, 0, 360,
                            "both",
                            "energyRadiusLog", "bar", configRFile, fileInfoFit, legendRL, 20, 1)
    legendRL.SetTextSize(0.045)
    legendRL.Draw()
    saveFigure(canvas, newPath, "energyMaxVsDistance" + str(configRFile.EnergyTime.Jump) + "Log")
    fileInfoFit.write("\n")

    fileInfoFit.write("---PEAK ENERGY VS DISTANCE \n")
    canvas = ROOT.TCanvas("energyRPeak")
    legendPR = ROOT.TLegend(0.75, 0.55, 0.9, 0.9)
    multiGraph = plotResult(selPixelsInfosPhi, 23, "distance [km]", configRFile.EnergyRadius.XMin,
                            configRFile.EnergyRadius.XMax,
                            configRFile.EnergyRadius.Jump, 13,
                            "energy peak [ADC]", configRFile.EnergyRadius.YMin, configRFile.EnergyRadius.YMax, 0, 360,
                            "both",
                            "energyRadius", "bar", configRFile, fileInfoFit, legendPR, 20, 1)
    legendPR.SetTextSize(0.04)
    legendPR.Draw()
    saveFigure(canvas, newPath, "energyPeakVsDistance" + str(configRFile.EnergyTime.Jump))
    fileInfoFit.write("\n")

    # in log scale:
    canvas = ROOT.TCanvas("energyPeakLog")
    gPad = ROOT.TPad("pad1", "pad1", 0.1, 0.1, 0.9, 0.9, 0)
    gPad.Draw()
    gPad.cd()
    gPad.SetLogy()
    legendPRL = ROOT.TLegend(0.7, 0.8, 0.9, 0.9)
    multiGraph = plotResult(selPixelsInfosPhi, 23, "distance [km]", configRFile.EnergyRadius.XMin,
                            configRFile.EnergyRadius.XMax,
                            configRFile.EnergyRadius.Jump, 13,
                            "energy peak [ADC]", configRFile.EnergyRadius.YMin, configRFile.EnergyRadius.YMax, 0, 360,
                            "both",
                            "energyRadiusLog", "bar", configRFile, fileInfoFit, legendPRL, 20, 1)
    legendPRL.SetTextSize(0.045)
    legendPRL.Draw()
    saveFigure(canvas, newPath, "energyPeakVsDistance" + str(configRFile.EnergyTime.Jump) + "Log")
    fileInfoFit.write("\n")
    fileInfoFit.close()


def plotEnergyResultsVsAngle(selPixelsInfosPhi, configRFile, newPath):

    if not os.path.exists(newPath):
        os.makedirs(newPath)
    fileInfoFit = open(newPath + "/infoFits.txt", "w+")

    # figure FrameAngles.png (histoPlot with blue and red sections designed)
    bestXCircleCenter = configRFile.CircleCoords.XCoord
    bestYCircleCenter = configRFile.CircleCoords.YCoord
    canvasFramePhi = ROOT.TCanvas()
    plotFramePhi = ROOT.TMultiGraph()
    frameToPlot = configRFile.EnergyAngle.FrameToPlot
    candidates = [i for i in range(len(selPixelsInfosPhi)) if ((selPixelsInfosPhi[i][6]-2 <= frameToPlot) and (selPixelsInfosPhi[i][6]+2 >= frameToPlot))]
    pointsToPlot = np.zeros((48, 48), dtype=np.float32)
    for t in range(len(candidates)):
            tt = candidates[t]
            pixX = selPixelsInfosPhi[tt][1]
            pixY = selPixelsInfosPhi[tt][2]
            pointsToPlot[int(pixX)][int(pixY)] = selPixelsInfosPhi[tt][13]
    (histo, canvasFramePhi) = fillHisto(newPath, configRFile.PixelFileName, configRFile.HH, pointsToPlot, 0, configRFile.EnergyAngle.MaxHisto, [bestXCircleCenter], [bestYCircleCenter], 1, "Frame %d " %(frameToPlot), "/FrameAngles.png")

    canvasCC = ROOT.TGraph()
    canvasCC.SetPoint(canvasCC.GetN(), bestXCircleCenter, bestYCircleCenter)
    canvasCC.SetPoint(canvasCC.GetN(), bestXCircleCenter + 400 * np.cos(
        configRFile.EnergyAngle.Angle1Min * math.pi / 180), bestYCircleCenter + 400 * np.sin(
        configRFile.EnergyAngle.Angle1Min * math.pi / 180))
    canvasCC.SetLineColor(2)
    canvasCC.SetLineStyle(1)
    canvasCC.SetLineWidth(4)
    canvasCCC = ROOT.TGraph()
    canvasCCC.SetPoint(canvasCCC.GetN(), bestXCircleCenter, bestYCircleCenter)
    canvasCCC.SetPoint(canvasCCC.GetN(), bestXCircleCenter + 400 * np.cos(
        configRFile.EnergyAngle.Angle1Max * math.pi / 180), bestYCircleCenter + 400 * np.sin(
        configRFile.EnergyAngle.Angle1Max * math.pi / 180))
    canvasCCC.SetLineColor(2)
    canvasCCC.SetLineStyle(1)
    canvasCCC.SetLineWidth(4)
    canvasCCCC = ROOT.TGraph()
    canvasCCCC.SetPoint(canvasCCCC.GetN(), bestXCircleCenter, bestYCircleCenter)
    canvasCCCC.SetPoint(canvasCCCC.GetN(), bestXCircleCenter + 400 * np.cos(
        configRFile.EnergyAngle.Angle2Min * math.pi / 180), bestYCircleCenter + 400 * np.sin(
        configRFile.EnergyAngle.Angle2Min * math.pi / 180))
    canvasCCCC.SetLineColor(4)
    canvasCCCC.SetLineStyle(1)
    canvasCCCC.SetLineWidth(4)
    canvasCCCCC = ROOT.TGraph()
    canvasCCCCC.SetPoint(canvasCCCCC.GetN(), bestXCircleCenter, bestYCircleCenter)
    canvasCCCCC.SetPoint(canvasCCCCC.GetN(), bestXCircleCenter + 400 * np.cos(
        configRFile.EnergyAngle.Angle2Max * math.pi / 180), bestYCircleCenter + 400 * np.sin(
        configRFile.EnergyAngle.Angle2Max * math.pi / 180))
    canvasCCCCC.SetLineColor(4)
    canvasCCCCC.SetLineStyle(1)
    canvasCCCCC.SetLineWidth(4)
    plotFramePhi.Add(canvasCC)
    plotFramePhi.Add(canvasCCC)
    plotFramePhi.Add(canvasCCCC)
    plotFramePhi.Add(canvasCCCCC)
    plotFramePhi.Draw("AL same")
    histo.Draw("COLZ same")
    canvasFramePhi.Update()
    plotFramePhi.GetXaxis().SetRangeUser(-150, 150)
    plotFramePhi.GetYaxis().SetRangeUser(-150, 150)
    saveFigure(canvasFramePhi, newPath, "FrameAngles")

    canvas = ROOT.TCanvas("energyMaxAngle")
    multiGraph = plotResult(selPixelsInfosPhi, 24, "angle [deg]", configRFile.EnergyAngle.XMin,
                            configRFile.EnergyAngle.XMax,
                            configRFile.EnergyAngle.Jump, 7,
                            "energy max [ADC]", configRFile.EnergyAngle.YMin, configRFile.EnergyAngle.YMax, 0, 360,
                            "both",
                            "none", "none", configRFile, fileInfoFit, ROOT.TLegend(), 20, 1)
    saveFigure(canvas, newPath, "energyMaxVsAngle" + str(configRFile.EnergyAngle.Jump))

    canvas = ROOT.TCanvas("energyMaxAngleVsTime")
    multiGraph = ROOT.TMultiGraph()
    legendAT = ROOT.TLegend(0.3, 0.8, 0.9, 0.9)
    multiGraph1 = plotResult(selPixelsInfosPhi, 6, "time [GTU]", configRFile.EnergyTime.XMin,
                             configRFile.EnergyTime.XMax,
                             configRFile.EnergyTime.Jump, 7,
                             "energy max [ADC]", configRFile.EnergyTime.YMin, configRFile.EnergyTime.YMax,
                             configRFile.EnergyAngle.Angle1Min, configRFile.EnergyAngle.Angle1Max,
                             "both",
                             "none", "none", configRFile, fileInfoFit, ROOT.TLegend(), 43, 2)
    multiGraph2 = plotResult(selPixelsInfosPhi, 6, "time [GTU]", configRFile.EnergyTime.XMin,
                             configRFile.EnergyTime.XMax,
                             configRFile.EnergyTime.Jump, 7,
                             "energy max [ADC]", configRFile.EnergyTime.YMin, configRFile.EnergyTime.YMax,
                             configRFile.EnergyAngle.Angle2Min, configRFile.EnergyAngle.Angle2Max,
                             "both",
                             "none", "none", configRFile, fileInfoFit, ROOT.TLegend(), 33, 4)
    multiGraph.Add(multiGraph1)
    multiGraph.Add(multiGraph2)
    multiGraph.Draw()
    tGraph = ROOT.TGraph()
    tGraph.SetMarkerStyle(20)
    tGraph.SetMarkerSize(0.5)
    tGraph.SetMarkerColor(2)
    legendAT.AddEntry(tGraph,
                      f"Sector between {configRFile.EnergyAngle.Angle1Min} and {configRFile.EnergyAngle.Angle1Max} degrees",
                      "p")
    tGraph2 = ROOT.TGraph()
    tGraph2.SetMarkerStyle(20)
    tGraph2.SetMarkerSize(0.5)
    tGraph2.SetMarkerColor(4)
    legendAT.AddEntry(tGraph2,
                      f"Sector between {configRFile.EnergyAngle.Angle2Min} and {configRFile.EnergyAngle.Angle2Max} degrees",
                      "p")
    legendAT.SetTextSize(0.04)
    legendAT.Draw()
    saveFigure(canvas, newPath, "energyMaxVsTime" + str(configRFile.EnergyAngle.Jump) + "Sectors")

    canvas = ROOT.TCanvas("energyMaxAngleVsTimeLog")
    gPad = ROOT.TPad("pad1", "pad1", 0.1, 0.1, 0.9, 0.9, 0)
    gPad.Draw()
    gPad.cd()
    gPad.SetLogy()
    multiGraph = ROOT.TMultiGraph()
    legendAT = ROOT.TLegend(0.3, 0.8, 0.9, 0.9)
    multiGraph1 = plotResult(selPixelsInfosPhi, 6, "time [GTU]", configRFile.EnergyTime.XMin,
                             configRFile.EnergyTime.XMax,
                             configRFile.EnergyTime.Jump, 7,
                             "energy max [ADC]", configRFile.EnergyTime.YMin, configRFile.EnergyTime.YMax,
                             configRFile.EnergyAngle.Angle1Min, configRFile.EnergyAngle.Angle1Max,
                             "both",
                             "none", "none", configRFile, fileInfoFit, ROOT.TLegend(), 43, 2)
    multiGraph2 = plotResult(selPixelsInfosPhi, 6, "time [GTU]", configRFile.EnergyTime.XMin,
                             configRFile.EnergyTime.XMax,
                             configRFile.EnergyTime.Jump, 7,
                             "energy max [ADC]", configRFile.EnergyTime.YMin, configRFile.EnergyTime.YMax,
                             configRFile.EnergyAngle.Angle2Min, configRFile.EnergyAngle.Angle2Max,
                             "both",
                             "none", "none", configRFile, fileInfoFit, ROOT.TLegend(), 33, 4)
    multiGraph.Add(multiGraph1)
    multiGraph.Add(multiGraph2)
    multiGraph.Draw()
    tGraph = ROOT.TGraph()
    tGraph.SetMarkerStyle(20)
    tGraph.SetMarkerSize(0.5)
    tGraph.SetMarkerColor(2)
    legendAT.AddEntry(tGraph,
                      f"Sector between {configRFile.EnergyAngle.Angle1Min} and {configRFile.EnergyAngle.Angle1Max} degrees",
                      "p")
    tGraph2 = ROOT.TGraph()
    tGraph2.SetMarkerStyle(20)
    tGraph2.SetMarkerSize(0.5)
    tGraph2.SetMarkerColor(4)
    legendAT.AddEntry(tGraph2,
                      f"Sector between {configRFile.EnergyAngle.Angle2Min} and {configRFile.EnergyAngle.Angle2Max} degrees",
                      "p")
    legendAT.SetTextSize(0.04)
    legendAT.Draw()
    saveFigure(canvas, newPath, "energyMaxVsTime" + str(configRFile.EnergyAngle.Jump) + "SectorsLog")

    canvas = ROOT.TCanvas("energyMaxAngleVsDistance")
    multiGraph = ROOT.TMultiGraph()
    legendAT = ROOT.TLegend(0.3, 0.8, 0.9, 0.9)
    multiGraph1 = plotResult(selPixelsInfosPhi, 23, "distance [km]", configRFile.EnergyRadius.XMin,
                             configRFile.EnergyRadius.XMax,
                             configRFile.EnergyRadius.Jump, 7,
                             "energy max [ADC]", configRFile.EnergyRadius.YMin, configRFile.EnergyRadius.YMax,
                             configRFile.EnergyAngle.Angle1Min, configRFile.EnergyAngle.Angle1Max,
                             "both",
                             "none", "none", configRFile, fileInfoFit, ROOT.TLegend(), 43, 2)
    multiGraph2 = plotResult(selPixelsInfosPhi, 23, "distance [km]", configRFile.EnergyRadius.XMin,
                             configRFile.EnergyRadius.XMax,
                             configRFile.EnergyRadius.Jump, 7,
                             "energy max [ADC]", configRFile.EnergyRadius.YMin, configRFile.EnergyRadius.YMax,
                             configRFile.EnergyAngle.Angle2Min, configRFile.EnergyAngle.Angle2Max,
                             "both",
                             "none", "none", configRFile, fileInfoFit, ROOT.TLegend(), 33, 4)
    multiGraph.Add(multiGraph1)
    multiGraph.Add(multiGraph2)
    multiGraph.Draw()
    tGraph = ROOT.TGraph()
    tGraph.SetMarkerStyle(20)
    tGraph.SetMarkerSize(0.5)
    tGraph.SetMarkerColor(2)
    legendAT.AddEntry(tGraph,
                      f"Sector between {configRFile.EnergyAngle.Angle1Min} and {configRFile.EnergyAngle.Angle1Max} degrees",
                      "p")
    tGraph2 = ROOT.TGraph()
    tGraph2.SetMarkerStyle(20)
    tGraph2.SetMarkerSize(0.5)
    tGraph2.SetMarkerColor(4)
    legendAT.AddEntry(tGraph2,
                      f"Sector between {configRFile.EnergyAngle.Angle2Min} and {configRFile.EnergyAngle.Angle2Max} degrees",
                      "p")
    legendAT.SetTextSize(0.04)
    legendAT.Draw()
    saveFigure(canvas, newPath, "energyMaxVsDistance" + str(configRFile.EnergyAngle.Jump) + "Sectors")

    canvas = ROOT.TCanvas("energyMaxAngleVsDistanceLog")
    gPad = ROOT.TPad("pad1", "pad1", 0.1, 0.1, 0.9, 0.9, 0)
    gPad.Draw()
    gPad.cd()
    gPad.SetLogy()
    multiGraph = ROOT.TMultiGraph()
    legendAT = ROOT.TLegend(0.3, 0.8, 0.9, 0.9)
    multiGraph1 = plotResult(selPixelsInfosPhi, 23, "distance [km]", configRFile.EnergyRadius.XMin,
                             configRFile.EnergyRadius.XMax,
                             configRFile.EnergyRadius.Jump, 7,
                             "energy max [ADC]", configRFile.EnergyRadius.YMin, configRFile.EnergyRadius.YMax,
                             configRFile.EnergyAngle.Angle1Min, configRFile.EnergyAngle.Angle1Max,
                             "both",
                             "none", "none", configRFile, fileInfoFit, ROOT.TLegend(), 43, 2)
    multiGraph2 = plotResult(selPixelsInfosPhi, 23, "distance [km]", configRFile.EnergyRadius.XMin,
                             configRFile.EnergyRadius.XMax,
                             configRFile.EnergyRadius.Jump, 7,
                             "energy max [ADC]", configRFile.EnergyRadius.YMin, configRFile.EnergyRadius.YMax,
                             configRFile.EnergyAngle.Angle2Min, configRFile.EnergyAngle.Angle2Max,
                             "both",
                             "none", "none", configRFile, fileInfoFit, ROOT.TLegend(), 33, 4)
    multiGraph.Add(multiGraph1)
    multiGraph.Add(multiGraph2)
    multiGraph.Draw()
    tGraph = ROOT.TGraph()
    tGraph.SetMarkerStyle(20)
    tGraph.SetMarkerSize(0.5)
    tGraph.SetMarkerColor(2)
    legendAT.AddEntry(tGraph,
                      f"Sector between {configRFile.EnergyAngle.Angle1Min} and {configRFile.EnergyAngle.Angle1Max} degrees",
                      "p")
    tGraph2 = ROOT.TGraph()
    tGraph2.SetMarkerStyle(20)
    tGraph2.SetMarkerSize(0.5)
    tGraph2.SetMarkerColor(4)
    legendAT.AddEntry(tGraph2,
                      f"Sector between {configRFile.EnergyAngle.Angle2Min} and {configRFile.EnergyAngle.Angle2Max} degrees",
                      "p")
    legendAT.SetTextSize(0.04)
    legendAT.Draw()
    saveFigure(canvas, newPath, "energyMaxVsDistance" + str(configRFile.EnergyAngle.Jump) + "SectorsLog")

    canvas = ROOT.TCanvas("energyPeakAngle")
    multiGraph = plotResult(selPixelsInfosPhi, 24, "angle [deg]", configRFile.EnergyAngle.XMin,
                            configRFile.EnergyAngle.XMax,
                            configRFile.EnergyAngle.Jump, 13,
                            "energy peak [ADC]", configRFile.EnergyAngle.YMin, configRFile.EnergyAngle.YMax, 0, 360,
                            "both",
                            "none", "none", configRFile, fileInfoFit, ROOT.TLegend(), 20, 1)
    canvas.SaveAs(newPath + "energyPeakVsAngle" + str(configRFile.EnergyAngle.Jump) + ".png")

    canvas = ROOT.TCanvas("energyPeakAngleVsTime")
    multiGraph = ROOT.TMultiGraph()
    legendAT = ROOT.TLegend(0.3, 0.8, 0.9, 0.9)
    multiGraph1 = plotResult(selPixelsInfosPhi, 6, "time [GTU]", configRFile.EnergyTime.XMin,
                             configRFile.EnergyTime.XMax,
                             configRFile.EnergyTime.Jump, 13,
                             "energy peak [ADC]", configRFile.EnergyTime.YMin, configRFile.EnergyTime.YMax,
                             configRFile.EnergyAngle.Angle1Min, configRFile.EnergyAngle.Angle1Max,
                             "both",
                             "none", "none", configRFile, fileInfoFit, ROOT.TLegend(), 43, 2)
    multiGraph2 = plotResult(selPixelsInfosPhi, 6, "time [GTU]", configRFile.EnergyTime.XMin,
                             configRFile.EnergyTime.XMax,
                             configRFile.EnergyTime.Jump, 13,
                             "energy peak [ADC]", configRFile.EnergyTime.YMin, configRFile.EnergyTime.YMax,
                             configRFile.EnergyAngle.Angle2Min, configRFile.EnergyAngle.Angle2Max,
                             "both",
                             "none", "none", configRFile, fileInfoFit, ROOT.TLegend(), 33, 4)
    multiGraph.Add(multiGraph1)
    multiGraph.Add(multiGraph2)
    multiGraph.Draw()
    tGraph = ROOT.TGraph()
    tGraph.SetMarkerStyle(20)
    tGraph.SetMarkerSize(0.5)
    tGraph.SetMarkerColor(2)
    legendAT.AddEntry(tGraph,
                      f"Sector between {configRFile.EnergyAngle.Angle1Min} and {configRFile.EnergyAngle.Angle1Max} degrees",
                      "p")
    tGraph2 = ROOT.TGraph()
    tGraph2.SetMarkerStyle(20)
    tGraph2.SetMarkerSize(0.5)
    tGraph2.SetMarkerColor(4)
    legendAT.AddEntry(tGraph2,
                      f"Sector between {configRFile.EnergyAngle.Angle2Min} and {configRFile.EnergyAngle.Angle2Max} degrees",
                      "p")
    legendAT.SetTextSize(0.04)
    legendAT.Draw()
    saveFigure(canvas, newPath, "energyPeakVsTime" + str(configRFile.EnergyAngle.Jump) + "Sectors")

    canvas = ROOT.TCanvas("energyPeakAngleVsTimeLog")
    gPad = ROOT.TPad("pad1", "pad1", 0.1, 0.1, 0.9, 0.9, 0)
    gPad.Draw()
    gPad.cd()
    gPad.SetLogy()
    multiGraph = ROOT.TMultiGraph()
    legendAT = ROOT.TLegend(0.3, 0.8, 0.9, 0.9)
    multiGraph1 = plotResult(selPixelsInfosPhi, 6, "time [GTU]", configRFile.EnergyTime.XMin,
                             configRFile.EnergyTime.XMax,
                             configRFile.EnergyTime.Jump, 13,
                             "energy peak [ADC]", configRFile.EnergyTime.YMin, configRFile.EnergyTime.YMax,
                             configRFile.EnergyAngle.Angle1Min, configRFile.EnergyAngle.Angle1Max,
                             "both",
                             "none", "none", configRFile, fileInfoFit, ROOT.TLegend(), 43, 2)
    multiGraph2 = plotResult(selPixelsInfosPhi, 6, "time [GTU]", configRFile.EnergyTime.XMin,
                             configRFile.EnergyTime.XMax,
                             configRFile.EnergyTime.Jump, 13,
                             "energy peak [ADC]", configRFile.EnergyTime.YMin, configRFile.EnergyTime.YMax,
                             configRFile.EnergyAngle.Angle2Min, configRFile.EnergyAngle.Angle2Max,
                             "both",
                             "none", "none", configRFile, fileInfoFit, ROOT.TLegend(), 33, 4)
    multiGraph.Add(multiGraph1)
    multiGraph.Add(multiGraph2)
    multiGraph.Draw()
    tGraph = ROOT.TGraph()
    tGraph.SetMarkerStyle(20)
    tGraph.SetMarkerSize(0.5)
    tGraph.SetMarkerColor(2)
    legendAT.AddEntry(tGraph,
                      f"Sector between {configRFile.EnergyAngle.Angle1Min} and {configRFile.EnergyAngle.Angle1Max} degrees",
                      "p")
    tGraph2 = ROOT.TGraph()
    tGraph2.SetMarkerStyle(20)
    tGraph2.SetMarkerSize(0.5)
    tGraph2.SetMarkerColor(4)
    legendAT.AddEntry(tGraph2,
                      f"Sector between {configRFile.EnergyAngle.Angle2Min} and {configRFile.EnergyAngle.Angle2Max} degrees",
                      "p")
    legendAT.SetTextSize(0.04)
    legendAT.Draw()
    saveFigure(canvas, newPath, "energyPeakVsTime" + str(configRFile.EnergyAngle.Jump) + "SectorsLog")

    canvas = ROOT.TCanvas("energyPeakAngleVsDistance")
    multiGraph = ROOT.TMultiGraph()
    legendAT = ROOT.TLegend(0.3, 0.8, 0.9, 0.9)
    multiGraph1 = plotResult(selPixelsInfosPhi, 23, "distance [km]", configRFile.EnergyRadius.XMin,
                             configRFile.EnergyRadius.XMax,
                             configRFile.EnergyRadius.Jump, 13,
                             "energy peak [ADC]", configRFile.EnergyRadius.YMin, configRFile.EnergyRadius.YMax,
                             configRFile.EnergyAngle.Angle1Min, configRFile.EnergyAngle.Angle1Max,
                             "both",
                             "none", "none", configRFile, fileInfoFit, ROOT.TLegend(), 43, 2)
    multiGraph2 = plotResult(selPixelsInfosPhi, 23, "distance [km]", configRFile.EnergyRadius.XMin,
                             configRFile.EnergyRadius.XMax,
                             configRFile.EnergyRadius.Jump, 13,
                             "energy peak [ADC]", configRFile.EnergyRadius.YMin, configRFile.EnergyRadius.YMax,
                             configRFile.EnergyAngle.Angle2Min, configRFile.EnergyAngle.Angle2Max,
                             "both",
                             "none", "none", configRFile, fileInfoFit, ROOT.TLegend(), 33, 4)
    multiGraph.Add(multiGraph1)
    multiGraph.Add(multiGraph2)
    multiGraph.Draw()
    tGraph = ROOT.TGraph()
    tGraph.SetMarkerStyle(20)
    tGraph.SetMarkerSize(0.5)
    tGraph.SetMarkerColor(2)
    legendAT.AddEntry(tGraph,
                      f"Sector between {configRFile.EnergyAngle.Angle1Min} and {configRFile.EnergyAngle.Angle1Max} degrees",
                      "p")
    tGraph2 = ROOT.TGraph()
    tGraph2.SetMarkerStyle(20)
    tGraph2.SetMarkerSize(0.5)
    tGraph2.SetMarkerColor(4)
    legendAT.AddEntry(tGraph2,
                      f"Sector between {configRFile.EnergyAngle.Angle2Min} and {configRFile.EnergyAngle.Angle2Max} degrees",
                      "p")
    legendAT.SetTextSize(0.04)
    legendAT.Draw()
    saveFigure(canvas, newPath, "energyPeakVsDistance" + str(configRFile.EnergyAngle.Jump) + "Sectors")

    canvas = ROOT.TCanvas("energyPeakAngleVsDistanceLog")
    gPad = ROOT.TPad("pad1", "pad1", 0.1, 0.1, 0.9, 0.9, 0)
    gPad.Draw()
    gPad.cd()
    gPad.SetLogy()
    multiGraph = ROOT.TMultiGraph()
    legendAT = ROOT.TLegend(0.3, 0.8, 0.9, 0.9)
    multiGraph1 = plotResult(selPixelsInfosPhi, 23, "distance [km]", configRFile.EnergyRadius.XMin,
                             configRFile.EnergyRadius.XMax,
                             configRFile.EnergyRadius.Jump, 13,
                             "energy peak [ADC]", configRFile.EnergyRadius.YMin, configRFile.EnergyRadius.YMax,
                             configRFile.EnergyAngle.Angle1Min, configRFile.EnergyAngle.Angle1Max,
                             "both",
                             "none", "none", configRFile, fileInfoFit, ROOT.TLegend(), 43, 2)
    multiGraph2 = plotResult(selPixelsInfosPhi, 23, "distance [km]", configRFile.EnergyRadius.XMin,
                             configRFile.EnergyRadius.XMax,
                             configRFile.EnergyRadius.Jump, 13,
                             "energy peak [ADC]", configRFile.EnergyRadius.YMin, configRFile.EnergyRadius.YMax,
                             configRFile.EnergyAngle.Angle2Min, configRFile.EnergyAngle.Angle2Max,
                             "both",
                             "none", "none", configRFile, fileInfoFit, ROOT.TLegend(), 33, 4)
    multiGraph.Add(multiGraph1)
    multiGraph.Add(multiGraph2)
    multiGraph.Draw()
    tGraph = ROOT.TGraph()
    tGraph.SetMarkerStyle(20)
    tGraph.SetMarkerSize(0.5)
    tGraph.SetMarkerColor(2)
    legendAT.AddEntry(tGraph,
                      f"Sector between {configRFile.EnergyAngle.Angle1Min} and {configRFile.EnergyAngle.Angle1Max} degrees",
                      "p")
    tGraph2 = ROOT.TGraph()
    tGraph2.SetMarkerStyle(20)
    tGraph2.SetMarkerSize(0.5)
    tGraph2.SetMarkerColor(4)
    legendAT.AddEntry(tGraph2,
                      f"Sector between {configRFile.EnergyAngle.Angle2Min} and {configRFile.EnergyAngle.Angle2Max} degrees",
                      "p")
    legendAT.SetTextSize(0.04)
    legendAT.Draw()
    saveFigure(canvas, newPath, "energyPeakVsDistance" + str(configRFile.EnergyAngle.Jump) + "SectorsLog")

    fileInfoFit.close()


def plotResults(pixelsInfosPhi, configRFile, newPath):
    if not os.path.exists(newPath):
        os.makedirs(newPath)

    ### SELECT RING:
    selPixelsInfosPhi = []
    for ll in range(len(pixelsInfosPhi)):
        if pixelsInfosPhi[ll][5] == 1:
            selPixelsInfosPhi.append(pixelsInfosPhi[ll])

    ### ENERGY:

    # energy vs time:
    plotEnergyResultsVsTime(selPixelsInfosPhi, configRFile, newPath + "vsTime/")

    # energy vs distance:
    plotEnergyResultsVsDistance(selPixelsInfosPhi, configRFile, newPath + "vsDistance/")

    # energy vs angle:
    plotEnergyResultsVsAngle(selPixelsInfosPhi, configRFile, newPath + "vsAngle/")



def plotResult(data, xLabel, xName, xLimBot, xLimTop, jumpX, yLabel, yName, yLimBot, yLimTop, angleBot, angleTop, plotType, fitType,
               fitPoints, configRFile, fileInfoFit, legend, colorMarkerScatter, colorMarkerBar):

    # sort data:
    dataSorted = np.zeros((len(data), len(data[0][:])))
    for l in range(len(data)):
        dataSorted[l][:] = np.array(data[l][:])
    dataSorted = dataSorted[dataSorted[:, xLabel].argsort()]
    data = dataSorted
    multiGraph = ROOT.TMultiGraph()

    # get arrays to plot:
    xArrayScatter = []
    yArrayScatter = []
    for ll in range(len(data)):
        if (xLimBot <= data[ll][xLabel] <= xLimTop) and (yLimBot <= data[ll][yLabel] <= yLimTop) and (angleBot <= data[ll][24] <= angleTop):
            xArrayScatter.append(data[ll][xLabel])
            yArrayScatter.append(data[ll][yLabel])
    result = ROOT.TGraph()
    for i in range(len(xArrayScatter)):
        result.SetPoint(result.GetN(), xArrayScatter[i], yArrayScatter[i])
    result.SetMarkerStyle(20)
    result.SetMarkerColor(colorMarkerScatter)
    result.SetMarkerSize(0.5)
    result.SetLineWidth(0)

    xArrayBar = []
    yArrayBar = []
    yErrArray = []
    zeros = []
    startCount = round(xLimBot / jumpX) * jumpX
    for ii in range(startCount, xLimTop + 1, jumpX):
        candidates = [i for i in range(len(data)) if
                      (data[i][xLabel] > ii - jumpX) and (data[i][xLabel] < ii + jumpX) and (
                                  yLimBot <= data[i][yLabel] <= yLimTop) and (angleBot <= data[i][24] <= angleTop)]
        if len(candidates) > 0:
            xArrayBar.append(ii)
            interval = np.array(candidates, dtype=int)
            yArrayBar.append(np.average(data[interval, yLabel]))
            yErrArray.append(np.std(data[interval, yLabel]))
            zeros.append(0)
    resultErr = ROOT.TGraphErrors(len(xArrayBar), array('f', xArrayBar), array('f', yArrayBar), array('f', zeros),
                                  array('f', yErrArray))
    resultErr.SetMarkerStyle(20)
    resultErr.SetMarkerColor(colorMarkerBar)
    resultErr.SetLineColor(1)
    resultErr.SetMarkerSize(0.5)

    # plot markers:
    if plotType == "scatter":
        multiGraph.Add(result, "AP")
    if plotType == "bar":
        multiGraph.Add(resultErr, "AP")
    if plotType == "both":
        multiGraph.Add(result, "AP")
        multiGraph.Add(resultErr, "AP")
    multiGraph.Draw("AP")

    # get fit line:
    if fitType == "radius":
        fileInfoFit.write(f"Time step: {configRFile.Radius.Jump} \n")
        # fit radius:
        configP = configRFile.Radius.RadiusSqrt
        function = RadiusFunction
        functionPlot = plotRadiusFunction
        functionR = getRSquaredRadius
        fileInfoFit.write("-- Results radius sqrt: \n")
        configPP = configP.Parameters
        initValues = [configPP.Velocity.Value, configPP.InitTime.Value, configPP.Altitude.Value]
        initBounds = ((configPP.Velocity.Min, configPP.InitTime.Min, configPP.Altitude.Min),
                      (configPP.Velocity.Max, configPP.InitTime.Max, configPP.Altitude.Max))
        fitLineRadius = fitResults(configP, function, functionR, functionPlot, initValues, initBounds, fileInfoFit,
                                   fitPoints, xArrayScatter, yArrayScatter, xArrayBar, yArrayBar, 2)
        multiGraph.Add(fitLineRadius, "AL")
        legend.AddEntry(fitLineRadius, "R(t)=#sqrt{v^{2}(t-t_{f})^{2}-(H-h_{0})^{2}}", "l")
        # fit linear radius:
        configP = configRFile.Radius.RadiusLinear
        function = RadiusFunctionLinear
        functionPlot = plotRadiusFunctionLinear
        functionR = getRSquaredRadiusLinear
        fileInfoFit.write("-- Results radius linear: \n")
        configPP = configP.Parameters
        initValues = [configPP.Velocity.Value, configPP.InitTime.Value]
        initBounds = ((configPP.Velocity.Min, configPP.InitTime.Min),
                      (configPP.Velocity.Max, configPP.InitTime.Max))
        fitLineRadiusLinear = fitResults(configP, function, functionR, functionPlot, initValues, initBounds,
                                         fileInfoFit,
                                         fitPoints, xArrayScatter, yArrayScatter, xArrayBar, yArrayBar, 3)
        multiGraph.Add(fitLineRadiusLinear, "AL")
        legend.AddEntry(fitLineRadiusLinear, "R(t)=v(t-t_{f})", "l")

    if fitType == "energyTime":
        fileInfoFit.write(f"Jump energy: {configRFile.EnergyTime.Jump} \n")
        # fit energy 1:
        configP = configRFile.EnergyTime.Energy2
        function = EnergyFunction
        functionPlot = plotEnergyFunction
        functionR = getRSquaredEnergy
        fileInfoFit.write("Results energy E(t) as r^(-2): \n")
        configPP = configP.Parameters
        initValues = [configPP.Velocity.Value, configPP.InitTime.Value, configPP.Altitude.Value, configPP.E1.Value]
        initBounds = ((configPP.Velocity.Min, configPP.InitTime.Min, configPP.Altitude.Min, configPP.E1.Min),
                      (configPP.Velocity.Max, configPP.InitTime.Max, configPP.Altitude.Max, configPP.E1.Max))
        fitLineEnergy1 = fitResults(configP, function, functionR, functionPlot, initValues, initBounds, fileInfoFit,
                                    fitPoints, xArrayScatter, yArrayScatter, xArrayBar, yArrayBar, 2)
        multiGraph.Add(fitLineEnergy1, "AL")
        legend.AddEntry(fitLineEnergy1, "E(t)=#frac{E_{1}}{v^{2}(t-t_{f})^{2}}", "l")
        # fit energy 2:
        configP = configRFile.EnergyTime.Energy32
        function = EnergyFunction2
        functionPlot = plotEnergyFunction2
        functionR = getRSquaredEnergy2
        fileInfoFit.write("Results energy E(t) as r^(-3/2): \n")
        configPP = configP.Parameters
        initValues = [configPP.Velocity.Value, configPP.InitTime.Value, configPP.Altitude.Value, configPP.E1.Value]
        initBounds = ((configPP.Velocity.Min, configPP.InitTime.Min, configPP.Altitude.Min, configPP.E1.Min),
                      (configPP.Velocity.Max, configPP.InitTime.Max, configPP.Altitude.Max, configPP.E1.Max))
        fitLineEnergy2 = fitResults(configP, function, functionR, functionPlot, initValues, initBounds, fileInfoFit,
                                    fitPoints, xArrayScatter, yArrayScatter, xArrayBar, yArrayBar, 3)
        multiGraph.Add(fitLineEnergy2, "AL")
        legend.AddEntry(fitLineEnergy2, "E(t)=#frac{E_{1}}{v^{3/2}(t-t_{f})^{3/2}}", "l")
        # fit energy 3:
        configP = configRFile.EnergyTime.EnergyLinear
        function = EnergyFunction3
        functionPlot = plotEnergyFunction3
        functionR = getRSquaredEnergy3
        fileInfoFit.write("Results energy E(t) as r: \n")
        configPP = configP.Parameters
        initValues = [configPP.Velocity.Value, configPP.InitTime.Value, configPP.Altitude.Value, configPP.E1.Value]
        initBounds = ((configPP.Velocity.Min, configPP.InitTime.Min, configPP.Altitude.Min, configPP.E1.Min),
                      (configPP.Velocity.Max, configPP.InitTime.Max, configPP.Altitude.Max, configPP.E1.Max))
        fitLineEnergy3 = fitResults(configP, function, functionR, functionPlot, initValues, initBounds, fileInfoFit,
                                    fitPoints, xArrayScatter, yArrayScatter, xArrayBar, yArrayBar, 4)
        multiGraph.Add(fitLineEnergy3, "AL")
        legend.AddEntry(fitLineEnergy3, "E(t)=E_{1}v(t-t_{f})", "l")

    if fitType == "energyRadius":
        fileInfoFit.write(f"Jump energy: {configRFile.EnergyRadius.Jump} \n")
        # fit energy 1:
        configP = configRFile.EnergyRadius.Energy2
        function = EnergyRFunction
        functionPlot = plotREnergyFunction
        functionR = getRSquaredEnergyR
        fileInfoFit.write("Results energy E(r) as r^(-2): \n")
        configPP = configP.Parameters
        initValues = [configPP.E1.Value]
        initBounds = ((configPP.E1.Min), (configPP.E1.Max))
        fitLineEnergyR1 = fitResults(configP, function, functionR, functionPlot, initValues, initBounds, fileInfoFit,
                                    fitPoints, xArrayScatter, yArrayScatter, xArrayBar, yArrayBar, 2)
        multiGraph.Add(fitLineEnergyR1, "AL")
        legend.AddEntry(fitLineEnergyR1, "E(r)=#frac{E_{1}}{(r)^{2}}", "l")
        # fit energy 2:
        configP = configRFile.EnergyRadius.Energy32
        function = EnergyRFunction2
        functionPlot = plotREnergyFunction2
        functionR = getRSquaredEnergyR2
        fileInfoFit.write("Results energy E(r) as r^(-3/2): \n")
        configPP = configP.Parameters
        initValues = [configPP.E1.Value]
        initBounds = ((configPP.E1.Min), (configPP.E1.Max))
        fitLineEnergyR2 = fitResults(configP, function, functionR, functionPlot, initValues, initBounds, fileInfoFit,
                                    fitPoints, xArrayScatter, yArrayScatter, xArrayBar, yArrayBar, 3)
        multiGraph.Add(fitLineEnergyR2, "AL")
        legend.AddEntry(fitLineEnergyR2, "E(t)=#frac{E_{1}}{(r)^{3/2}}", "l")
        # fit energy 3:
        configP = configRFile.EnergyRadius.EnergyLinear
        function = EnergyRFunction3
        functionPlot = plotREnergyFunction3
        functionR = getRSquaredEnergyR3
        fileInfoFit.write("Results energy E(r) as r: \n")
        configPP = configP.Parameters
        initValues = [configPP.E1.Value]
        initBounds = ((configPP.E1.Min), (configPP.E1.Max))
        fitLineEnergyR3 = fitResults(configP, function, functionR, functionPlot, initValues, initBounds, fileInfoFit,
                                    fitPoints, xArrayScatter, yArrayScatter, xArrayBar, yArrayBar, 4)
        multiGraph.Add(fitLineEnergyR3, "AL")
        legend.AddEntry(fitLineEnergyR3, "E(t)=E_{1}r", "l")

    if (fitType == "energyTimeLog") or (fitType == "energyRadiusLog"):
        # fit log:
        if fitType == "energyTimeLog":
            fileInfoFit.write("Results energy time log: \n")
            configP = configRFile.EnergyTime.EnergyLog
        else:
            fileInfoFit.write("Results energy distance log: \n")
            configP = configRFile.EnergyRadius.EnergyLog
        function = EnergyFunctionLog
        functionPlot = plotEnergyFunctionLog
        functionR = getRSquaredEnergyLog
        initValues = [0, 0]
        initBounds = ((-999, -999), (999, 999))
        fitLineEnergyLog = fitResults(configP, function, functionR, functionPlot, initValues, initBounds, fileInfoFit,
                                    fitPoints, xArrayScatter, yArrayScatter, xArrayBar, yArrayBar, 2)
        multiGraph.Add(fitLineEnergyLog, "AL")
        legend.AddEntry(fitLineEnergyLog, "E(x)=e^{ax+b}", "l")

    multiGraph.Draw()

    # save .png:
    multiGraph.GetXaxis().SetTitle(xName)
    multiGraph.GetYaxis().SetTitle(yName)
    multiGraph.GetXaxis().SetRangeUser(xLimBot - 10, xLimTop + 10)
    multiGraph.GetYaxis().SetRangeUser(1, yLimTop + 50)

    return multiGraph


def fitResults(configP, function, functionR, functionPlot, initValues, initBounds, fileInfoFit, fitPoints,
               xArrayScatter, yArrayScatter, xArrayBar, yArrayBar, color):
    pointsToFit = []
    arrayFitX = []
    arrayFitY = []
    if fitPoints == "scatter":
        pointsToFit = [i for i in range(len(xArrayScatter)) if (configP.XMin <= xArrayScatter[i] <= configP.XMax) and (
                    configP.YMin <= yArrayScatter[i] <= configP.YMax)]
        if len(pointsToFit) > 0:
            arrayFitX = np.array(xArrayScatter)[np.array(pointsToFit, dtype=int)]
            arrayFitY = np.array(yArrayScatter)[np.array(pointsToFit, dtype=int)]
    if fitPoints == "bar":
        pointsToFit = [i for i in range(len(xArrayBar)) if (configP.XMin <= xArrayBar[i] <= configP.XMax) and (
                    configP.YMin <= yArrayBar[i] <= configP.YMax)]
        if len(pointsToFit) > 0:
            arrayFitX = np.array(xArrayBar)[np.array(pointsToFit, dtype=int)]
            arrayFitY = np.array(yArrayBar)[np.array(pointsToFit, dtype=int)]

    if len(arrayFitX) > 0:

        popt, pcov = sciopt.curve_fit(function, arrayFitX, arrayFitY, initValues, bounds=initBounds)
        perr = np.sqrt(np.diag(pcov))
        rSquared = functionR(arrayFitX, arrayFitY, popt)

        fileInfoFit.write("R2: %5.4f \n" %(rSquared))
        fileInfoFit.write(f"Fit results: {popt} \n")
        fileInfoFit.write(f"Fit errors: {perr} \n")

        # plot fit line:
        plotFitLine = ROOT.TGraph()
        for r in range(configP.XMin, configP.XMax):
            plotFitLine.SetPoint(plotFitLine.GetN(), r, functionPlot(r, popt))
        plotFitLine.SetLineStyle(1)
        plotFitLine.SetLineColor(color)
        plotFitLine.SetLineWidth(1)

    return plotFitLine


def RadiusFunction(t, v, t0, h0):
    Hiono = 90
    tf = t0 - ((Hiono - h0) / v)
    t1 = (t - tf)
    return np.sqrt(np.power(v * t1, 2) - np.power(Hiono - h0, 2))


def RadiusFunctionLinear(t, v, t0):
    t1 = (t - t0)
    return v * t1


def EnergyFunction(t, v, t0, h0, E1):
    Hiono = 90
    tf = t0 - ((Hiono - h0) / v)
    t1 = (t - tf)
    return 1e06 * E1 / np.power(v * t1, 2)


def EnergyFunctionLog(x, E1, E2):
    return np.exp(E2) * np.exp(x * E1)


def EnergyFunction2(t, v, t0, h0, E1):
    Hiono = 90
    tf = t0 - ((Hiono - h0) / v)
    t1 = (t - tf)
    return np.power(10, 4.5) * E1 / np.power(v * t1, 3 / 2)


def EnergyFunction3(t, v, t0, h0, E1):
    Hiono = 90
    tf = t0 - ((Hiono - h0) / v)
    t1 = (t - tf)
    return E1 * (v * t1)


def EnergyRFunction(r, E1):
    return 1e06 * E1 / np.power(r, 2)


def EnergyRFunction2(r, E1):
    return np.power(10, 4.5) * E1 / np.power(r, 3 / 2)


def EnergyRFunction3(r, E1):
    return E1 * r


def plotRadiusFunction(t, popt):
    v = popt[0]
    t0 = popt[1]
    h0 = popt[2]
    Hiono = 90
    tf = t0 - ((Hiono - h0) / v)
    t1 = (t - tf)
    return np.sqrt(np.power(v * t1, 2) - np.power(Hiono - h0, 2))


def plotRadiusFunctionLinear(t, popt):
    v = popt[0]
    t0 = popt[1]
    t1 = (t - t0)
    return v * t1


def plotEnergyFunction(t, popt):
    v = popt[0]
    t0 = popt[1]
    h0 = popt[2]
    E1 = popt[3]
    Hiono = 90
    tf = t0 - ((Hiono - h0) / v)
    t1 = (t - tf)
    return 1e06 * E1 / np.power(v * t1, 2)


def plotEnergyFunctionLog(x, popt):
    E1 = popt[0]
    E2 = popt[1]
    return np.exp(E1*x) * np.exp(E2)


def plotEnergyFunction2(t, popt):
    v = popt[0]
    t0 = popt[1]
    h0 = popt[2]
    E1 = popt[3]
    Hiono = 90
    tf = t0 - ((Hiono - h0) / v)
    t1 = (t - tf)
    return np.power(10, 4.5) * E1 / np.power(v * t1, 3 / 2)


def plotEnergyFunction3(t, popt):
    v = popt[0]
    t0 = popt[1]
    h0 = popt[2]
    E1 = popt[3]
    Hiono = 90
    tf = t0 - ((Hiono - h0) / v)
    t1 = (t - tf)
    return E1 * (v * t1)


def plotREnergyFunction(r, popt):
    E1 = popt[0]
    return 1e06 * E1 / np.power(r, 2)


def plotREnergyFunction2(r, popt):
    E1 = popt[0]
    return np.power(10, 4.5) * E1 / np.power(r, (3 / 2))


def plotREnergyFunction3(r, popt):
    E1 = popt[0]
    return E1 * r


def getRSquaredRadius(arrayTime, arrayRadius, popt):
    meanRadius = np.average(arrayRadius)
    r_squared = 0
    ss_res = 0
    ss_tot = 0
    for l in range(len(arrayTime)):
        res = arrayRadius[l] - RadiusFunction(arrayTime[l], popt[0], popt[1], popt[2])
        ss_res = ss_res + res ** 2
        ss_tot = ss_tot + (arrayRadius[l] - meanRadius) ** 2
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared


def getRSquaredRadiusLinear(arrayTime, arrayRadius, popt):
    meanRadius = np.average(arrayRadius)
    r_squared = 0
    ss_res = 0
    ss_tot = 0
    for l in range(len(arrayTime)):
        res = arrayRadius[l] - RadiusFunctionLinear(arrayTime[l], popt[0], popt[1])
        ss_res = ss_res + res ** 2
        ss_tot = ss_tot + (arrayRadius[l] - meanRadius) ** 2
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared


def getRSquaredEnergy(arrayTime, arrayEnergy, popt):
    meanEnergy = np.average(arrayEnergy)
    r_squared = 0
    ss_res = 0
    ss_tot = 0
    for l in range(len(arrayTime)):
        res = arrayEnergy[l] - EnergyFunction(arrayTime[l], popt[0], popt[1], popt[2], popt[3])
        ss_res = ss_res + res ** 2
        ss_tot = ss_tot + (arrayEnergy[l] - meanEnergy) ** 2
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared


def getRSquaredEnergyLog(arrayX, arrayEnergy, popt):
    meanEnergy = np.average(arrayEnergy)
    r_squared = 0
    ss_res = 0
    ss_tot = 0
    for l in range(len(arrayX)):
        res = arrayEnergy[l] - EnergyFunctionLog(arrayX[l], popt[0], popt[1])
        ss_res = ss_res + res ** 2
        ss_tot = ss_tot + (arrayEnergy[l] - meanEnergy) ** 2
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared


def getRSquaredEnergy2(arrayTime, arrayEnergy, popt):
    meanEnergy = np.average(arrayEnergy)
    r_squared = 0
    ss_res = 0
    ss_tot = 0
    for l in range(len(arrayTime)):
        res = arrayEnergy[l] - EnergyFunction2(arrayTime[l], popt[0], popt[1], popt[2], popt[3])
        ss_res = ss_res + res ** 2
        ss_tot = ss_tot + (arrayEnergy[l] - meanEnergy) ** 2
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared


def getRSquaredEnergy3(arrayTime, arrayEnergy, popt):
    meanEnergy = np.average(arrayEnergy)
    r_squared = 0
    ss_res = 0
    ss_tot = 0
    for l in range(len(arrayTime)):
        res = arrayEnergy[l] - EnergyFunction3(arrayTime[l], popt[0], popt[1], popt[2], popt[3])
        ss_res = ss_res + res ** 2
        ss_tot = ss_tot + (arrayEnergy[l] - meanEnergy) ** 2
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared


def getRSquaredEnergyR(arrayDistance, arrayEnergy, popt):
    meanEnergy = np.average(arrayEnergy)
    r_squared = 0
    ss_res = 0
    ss_tot = 0
    for l in range(len(arrayDistance)):
        res = arrayEnergy[l] - EnergyRFunction(arrayDistance[l], popt[0])
        ss_res = ss_res + res ** 2
        ss_tot = ss_tot + (arrayEnergy[l] - meanEnergy) ** 2
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared


def getRSquaredEnergyR2(arrayDistance, arrayEnergy, popt):
    meanEnergy = np.average(arrayEnergy)
    r_squared = 0
    ss_res = 0
    ss_tot = 0
    for l in range(len(arrayDistance)):
        res = arrayEnergy[l] - EnergyRFunction2(arrayDistance[l], popt[0])
        ss_res = ss_res + res ** 2
        ss_tot = ss_tot + (arrayEnergy[l] - meanEnergy) ** 2
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared


def getRSquaredEnergyR3(arrayDistance, arrayEnergy, popt):
    meanEnergy = np.average(arrayEnergy)
    r_squared = 0
    ss_res = 0
    ss_tot = 0
    for l in range(len(arrayDistance)):
        res = arrayEnergy[l] - EnergyRFunction3(arrayDistance[l], popt[0])
        ss_res = ss_res + res ** 2
        ss_tot = ss_tot + (arrayEnergy[l] - meanEnergy) ** 2
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared
