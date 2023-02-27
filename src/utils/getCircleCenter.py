import os
import ROOT
import numpy as np
from src.utils.fillHisto import fillHisto


def getCircleCenter(newPath, pixelsInfos, data, configFile, frameStart):
    if not os.path.exists(newPath):
        os.makedirs(newPath)
    # open .txt file:
    fileDebug = open(newPath + "/infoCenter.txt", "w+")  # open debugFile
    fileDebugCoord = open(newPath + "/centerCoord.txt", "w+")  # open debugFile

    arrayTime = []
    xCentres = []
    yCentres = []
    arrayBestTime = []
    bestXCentres = []
    bestYCentres = []

    allCandidates = [i for i in range(len(pixelsInfos)) if pixelsInfos[i][5] == 1]

    # #if second peak:
    # secondPeakCandidates = []
    # for t in range(len(allCandidates)):
    #     tt = allCandidates[t]
    #     timeMax = pixelsInfos[tt][6]
    #     peaks = [x for x in range(len(pixelsInfos)) if (pixelsInfos[x][0] == pixelsInfos[tt][0]) and (pixelsInfos[x][5] > 1) and (pixelsInfos[x][6]-timeMax > 0)]
    #     if len(peaks) > 0:
    #         diffs = []
    #         for p in range(len(peaks)):
    #             diffs.append(abs(timeMax-pixelsInfos[peaks[p]][6]))
    #         minDiffs=min(diffs)
    #         secondP = [x for x in range(len(peaks)) if diffs[x] == minDiffs][0]
    #         secondPeakCandidates.append(peaks[secondP])
    # allCandidates = secondPeakCandidates

    # get array of centers:
    for time in range(frameStart, configFile.FrameEnd):

        candidates = [allCandidates[i] for i in range(len(allCandidates)) if (
                (pixelsInfos[allCandidates[i]][6] - 2 <= time) and (pixelsInfos[allCandidates[i]][6] + 2 >= time))]

        # plot pixels:
        plotPIXELS = ROOT.TCanvas("plotPixels" + str(time))
        plotPIX = ROOT.TMultiGraph()
        for t in range(len(candidates)):
            plotEnePixel = ROOT.TGraph()
            tt = candidates[t]
            print("%d \t %d \n" % (time, pixelsInfos[tt][0]))
            px = pixelsInfos[tt][1]
            py = pixelsInfos[tt][2]
            for timeP in range(len(data[0, 0, :, px, py])):
                plotEnePixel.SetPoint(plotEnePixel.GetN(), timeP + frameStart, data[0, 0, timeP, px, py])
            plotEnePixel.SetMarkerStyle(20)
            plotEnePixel.SetMarkerColor(t + 1)
            if t >= 9:
                plotEnePixel.SetMarkerColor(t + 2)
            plotEnePixel.SetMarkerSize(0.5)
            plotPIX.Add(plotEnePixel)
        plotPIX.Draw("APL")
        plotPIX.GetXaxis().SetTitle("Time")
        plotPIX.GetYaxis().SetTitle("Energy [ADC]")
        plotPIXELS.SaveAs(newPath + "/Pixels" + str(time) + ".png")
        if len(candidates) < 5:
            continue
        xCoord = np.array([pixelsInfos[p][3] for p in candidates])
        yCoord = np.array([pixelsInfos[p][4] for p in candidates])
        energies = np.array([pixelsInfos[p][13] for p in candidates])  # above half max energy
        (singleCentreX, singleCentreY, errCentreX, errCentreY, singleRadius) = getCircleCenterSingleFrame(xCoord,
                                                                                                          yCoord,
                                                                                                          energies,
                                                                                                          configFile.NominalCentreX,
                                                                                                          configFile.NominalCentreY,
                                                                                                          500)
        arrayTime.append(time)
        xCentres.append(singleCentreX)
        yCentres.append(singleCentreY)
        # plot single frame (energy):
        pointsToPlot = np.zeros((48, 48), dtype=np.float32)
        markX = []
        markY = []
        for t in range(len(candidates)):
            tt = candidates[t]
            pixX = pixelsInfos[tt][1]
            pixY = pixelsInfos[tt][2]
            pointsToPlot[pixX][pixY] = pixelsInfos[tt][13]  # above half max energy
        canvasH = ROOT.TCanvas("%d" % time)
        histoH = ROOT.TH2Poly()
        (histoH, canvasH) = fillHisto(newPath, configFile.PixelFileName, configFile.HH, pointsToPlot, 0,
                                      configFile.MaxHisto, [], [], 1, "",
                                      "/plotSingleFrame" + str(time) + ".png")
        histoH.Draw("COLZ")
        canvasH.Update()
        if (errCentreX < 10 and errCentreY < 10 and abs(singleCentreX) < 5000 and abs(
                singleCentreY) < 5000 and configFile.GoodFrameStart <= time <= configFile.GoodFrameEnd):
            bestXCentres.append(singleCentreX)
            bestYCentres.append(singleCentreY)
            arrayBestTime.append(time)
        fileDebugCoord.write("%d \t %5.2f \t %5.2f \t %5.2f \n" % (0, time, singleCentreX, singleCentreY))
        arc = ROOT.TArc(singleCentreX, singleCentreY, singleRadius)
        arc.SetLineColor(2)
        arc.SetLineWidth(2)
        arc.SetFillStyle(0)
        arc.Draw()
        mC = ROOT.TMarker(singleCentreX, singleCentreY, 20)
        mC.SetMarkerColor(2)
        mC.Draw()
        for k in range(len(markX)):
            m = ROOT.TMarker()
            m.SetMarkerStyle(20)
            m.SetMarkerColor(1)
            m.DrawMarker(markX[k], markY[k])
        fileDebugCoord.write("%d \t %5.2f \t %5.2f \t %5.2f \n" % (10, time, singleCentreX, singleCentreY))
        canvasH.SaveAs(newPath + "/plotSingleFrame" + str(time) + ".png")

    bestCentreXMean = 0
    bestCentreYMean = 0
    errCentreX = 0
    errCentreY = 0

    # get mean values (excluding outliers!):
    if len(bestXCentres) > 0:
        bestCentreXMean = np.mean(bestXCentres)
        bestCentreXStd = np.std(bestXCentres)
        noOutliersX = [bestXCentres[i] for i in range(len(bestXCentres)) if
                       abs(bestXCentres[i] - bestCentreXMean) < bestCentreXStd]
        bestCentreXMean = np.mean(noOutliersX)
        errCentreX = np.std(noOutliersX)
        bestCentreYMean = np.mean(bestYCentres)
        bestCentreYStd = np.std(bestYCentres)
        noOutliersY = [bestYCentres[i] for i in range(len(bestYCentres)) if
                       abs(bestYCentres[i] - bestCentreYMean) < bestCentreYStd]
        bestCentreYMean = np.mean(noOutliersY)
        errCentreY = np.std(noOutliersY)
    fileDebug.write("Centre Circle: \n X: %5.2f km \n ErrX: %3.2f km \n Y: %5.2f km \n ErrY: %3.2f km \n" % (
        bestCentreXMean, errCentreX, bestCentreYMean, errCentreY))
    fileDebug.write("Radial distance: %5.2f km \n" % (getDistance(bestCentreXMean, bestCentreYMean, 0, 0)))

    # plot result:
    plotPosX = ROOT.TGraph()
    for x in range(len(xCentres)):
        plotPosX.SetPoint(plotPosX.GetN(), arrayTime[x], xCentres[x])
    plotPosX.SetMarkerStyle(20)
    plotPosX.SetMarkerColor(4)
    plotPosX.SetMarkerSize(0.5)
    plotBestPosX = ROOT.TGraph()
    for x in range(len(bestXCentres)):
        plotBestPosX.SetPoint(plotBestPosX.GetN(), arrayBestTime[x], bestXCentres[x])
    plotBestPosX.SetMarkerStyle(20)
    plotBestPosX.SetMarkerColor(2)
    plotBestPosX.SetMarkerSize(0.5)
    plotPosY = ROOT.TGraph()
    for y in range(len(yCentres)):
        plotPosY.SetPoint(plotPosY.GetN(), arrayTime[y], yCentres[y])
    plotPosY.SetMarkerStyle(20)
    plotPosY.SetMarkerColor(4)
    plotPosY.SetMarkerSize(0.5)
    plotBestPosY = ROOT.TGraph()
    plotBestPosY = ROOT.TGraph()
    for y in range(len(bestYCentres)):
        plotBestPosY.SetPoint(plotBestPosY.GetN(), arrayBestTime[y], bestYCentres[y])
    plotBestPosY.SetMarkerStyle(20)
    plotBestPosY.SetMarkerColor(2)
    plotBestPosY.SetMarkerSize(0.5)
    plotPosXFitCircles = ROOT.TCanvas("plotPosXFitCircles")
    plotX = ROOT.TMultiGraph()
    plotX.Add(plotPosX)
    plotX.Add(plotBestPosX)
    plotX.Draw("AP")
    lineBest = ROOT.TLine(frameStart - 10, bestCentreXMean, frameStart + 500, bestCentreXMean)
    lineBest.SetLineColor(2)
    lineBest.SetLineStyle(9)
    lineBest.Draw()
    plotX.GetXaxis().SetTitle("time [#frame]")
    plotX.GetYaxis().SetTitle("X coordinate [km]")
    plotX.GetYaxis().SetRangeUser(bestCentreXMean - 500, bestCentreXMean + 500)
    plotPosXFitCircles.SaveAs(newPath + "/" + "plotPosXFitCircles.png")
    plotPosYFitCircles = ROOT.TCanvas("plotPosYFitCircles")
    plotY = ROOT.TMultiGraph()
    plotY.Add(plotPosY)
    plotY.Add(plotBestPosY)
    plotY.Draw("AP")
    lineBest = ROOT.TLine(frameStart - 10, bestCentreYMean, frameStart + 500, bestCentreYMean)
    lineBest.SetLineColor(2)
    lineBest.SetLineStyle(9)
    lineBest.Draw()
    plotY.GetXaxis().SetTitle("time [#frame]")
    plotY.GetYaxis().SetTitle("Y coordinate [km]")
    plotY.GetYaxis().SetRangeUser(bestCentreYMean - 500, bestCentreYMean + 500)
    plotPosYFitCircles.SaveAs(newPath + "/" + "plotPosYFitCircles.png")
    fileDebug.close()
    fileDebugCoord.close()
    return bestCentreXMean, bestCentreYMean, errCentreX, errCentreY


def getDistance(x0, y0, x1, y1):
    # get distance between (x0,y0) and (x1,y1);

    distance = 0
    distX = np.power(x1 - x0, 2)
    distY = np.power(y1 - y0, 2)
    distance = np.sqrt(distX + distY)
    return distance


def getCircleCenterSingleFrame(xCoord, yCoord, energies, nominalCentreX, nominalCentreY, maxSteps):
    centreX = 0
    centreY = 0
    errCentreX = 0
    errCentreY = 0
    energies = energies / sum(energies)  # normalize energies
    count = 0
    iteration = 0
    goodStep = -1
    arrayCentreX = []
    arrayCentreY = []
    lamda = 0.5
    lamdas = []
    lamdas.append(lamda)
    distances = getDistances(xCoord, yCoord, nominalCentreX, nominalCentreY)
    radius = np.median(distances)
    E = getErrorFunction(distances, radius)
    Errors = []
    Errors.append(E)
    # iteration:
    while count < maxSteps:
        iteration += 1
        # get new centre coordinates:
        (d1, d2) = getParametersD1D2(xCoord, yCoord, nominalCentreX, nominalCentreY, energies)
        (newCentreX, newCentreY) = updateCentreCoordinates(d1, d2, lamda, nominalCentreX, nominalCentreY)
        # get new values:
        newDistances = getDistances(xCoord, yCoord, newCentreX, newCentreY)
        newRadius = np.median(newDistances)
        newE = getErrorFunction(newDistances, newRadius)
        newLamda = updateLamda(E, newE, lamda)
        # update old values:
        arrayCentreX.append(newCentreX)
        arrayCentreY.append(newCentreY)
        lamdas.append(newLamda)
        Errors.append(newE)
        nominalCentreX = newCentreX
        nominalCentreY = newCentreY
        lamda = newLamda
        distances = newDistances
        E = newE
        if lamda >= 0.01:
            goodStep = iteration
            count = 0
        if (lamda < 0.01) & (goodStep > 0):  # first time lamda is zero
            count += 1
        if iteration >= 3000:  # force stop if not convergent
            count = maxSteps + 1
            goodStep = 1500
        if 0 != 1:
            pass
        else:  # plotting results
            plotXCentres = ROOT.TCanvas("plotXCentre")
            plotXCentre = ROOT.TGraph()
            for xx in range(1, len(arrayCentreX)):
                plotXCentre.SetPoint(plotXCentre.GetN(), xx, arrayCentreX[xx])
                plotXCentre.Draw()
            line = ROOT.TLine(goodStep, min(arrayCentreX) - 100, goodStep, max(arrayCentreX) + 100);
            line.SetLineColor(2)
            line.SetLineStyle(9)
            line.Draw()
            plotXCentre.SetTitle("X coordinate")
            plotXCentre.GetXaxis().SetTitle("Iteration")
            plotXCentre.GetYaxis().SetTitle("X coordinate [km]")
            plotXCentres.SaveAs("XCoordinateIteration.png")
            plotYCentres = ROOT.TCanvas("plotYCentre")
            plotYCentre = ROOT.TGraph()
            for yy in range(1, len(arrayCentreY)):
                plotYCentre.SetPoint(plotYCentre.GetN(), yy, arrayCentreY[yy])
                plotYCentre.Draw()
            line = ROOT.TLine(goodStep, min(arrayCentreY) - 100, goodStep, max(arrayCentreY) + 100);
            line.SetLineColor(2)
            line.SetLineStyle(9)
            line.Draw()
            plotYCentre.SetTitle("Y coordinate")
            plotYCentre.GetXaxis().SetTitle("Iteration")
            plotYCentre.GetYaxis().SetTitle("Y coordinate [km]")
            plotYCentres.SaveAs("YCoordinateIteration.png")
            plotLamdas = ROOT.TCanvas("plotLamda")
            plotLamda = ROOT.TGraph()
            for l in range(1, len(lamdas)):
                plotLamda.SetPoint(plotLamda.GetN(), l, lamdas[l])
                plotLamda.Draw()
            line = ROOT.TLine(goodStep, min(lamdas) - 100, goodStep, max(lamdas) + 100);
            line.SetLineColor(2)
            line.SetLineStyle(9)
            line.Draw()
            plotLamda.SetTitle("Lamda")
            plotLamda.GetXaxis().SetTitle("Iteration")
            plotLamda.GetYaxis().SetTitle("Lamda")
            plotLamdas.SaveAs("Lamda.png")
            plotErrorFunctions = ROOT.TCanvas("plotErrorFunction")
            plotErrorFunction = ROOT.TGraph()
            for e in range(1, len(Errors)):
                plotErrorFunction.SetPoint(plotErrorFunction.GetN(), e, Errors[e])
                plotErrorFunction.Draw()
            line = ROOT.TLine(goodStep, min(Errors) - 100, goodStep, max(Errors) + 100);
            line.SetLineColor(2)
            line.SetLineStyle(9)
            line.Draw()
            plotErrorFunction.SetTitle("Error function")
            plotErrorFunction.GetXaxis().SetTitle("Iteration")
            plotErrorFunction.GetYaxis().SetTitle("Error function")
            plotErrorFunctions.SaveAs("ErrorFunction.png")
        if goodStep > 0:
            centreX = np.mean(arrayCentreX[goodStep:])
            errCentreX = np.std(arrayCentreX[goodStep:])
            centreY = np.mean(arrayCentreY[goodStep:])
            errCentreY = np.std(arrayCentreY[goodStep:])
            finalDistances = getDistances(xCoord, yCoord, centreX, centreY)
            finalRadius = np.median(finalDistances)
    return centreX, centreY, errCentreX, errCentreY, finalRadius


def getDistances(x, y, centreX, centreY):
    distances = []
    if len(x) != len(y):
        print("Error in getting distances: array of coordinates don't match!")
        return distances
    for i in range(len(x)):
        distances.append(getDistance(x[i], y[i], centreX, centreY))
    return distances


def getErrorFunction(residuals, radius):
    # return the error function E;
    E = 0
    for i in range(len(residuals)):
        E += np.abs(residuals[i] - radius)
    return E


def getParametersD1D2(x, y, centreX, centreY, weights):
    d1 = 0
    d2 = 0
    cos = getCos(x, y, centreX, centreY, weights)
    sin = getSin(x, y, centreX, centreY, weights)
    distances = getDistances(x, y, centreX, centreY)
    radius = np.median(distances)
    (EaPlus, EbPlus, EaMinus, EbMinus) = getErrors(cos, sin, distances, radius)
    (alpha, beta) = getAlphaBeta(EaPlus, EbPlus, EaMinus, EbMinus)
    d1 = -(alpha * EaMinus + (1 - alpha) * EaPlus)
    d2 = -(beta * EbMinus + (1 - beta) * EbPlus)
    return d1, d2


def updateCentreCoordinates(d1, d2, lamda, centreX, centreY):
    # return the updated coordinates of the circle center;
    newCentreX = centreX + lamda * d1
    newCentreY = centreY + lamda * d2
    return newCentreX, newCentreY


def updateLamda(oldE, newE, oldLamda):
    # update lamda;
    if newE < oldE:
        newLamda = oldLamda * 1.1
    else:
        newLamda = oldLamda * 0.9
    return newLamda


def getCos(x, y, centreX, centreY, weights):
    # return the cos of the distances of the data points from the circle center;
    cos = []
    if len(x) != len(y):
        print("Error in getting cos: array of coordinates don't match!")
        return cos
    for i in range(len(x)):
        distX = weights[i] * (x[i] - centreX)
        distance = getDistance(x[i], y[i], centreX, centreY)
        cos.append(distX / distance)
    return cos


def getSin(x, y, centreX, centreY, weights):
    # return the sin of the distances of the data points from the circle center;
    sin = []
    if len(x) != len(y):
        print("Error in getting sin: array of coordinates don't match!")
        return sin
    for i in range(len(x)):
        distY = weights[i] * (y[i] - centreY)
        distance = getDistance(x[i], y[i], centreX, centreY)
        sin.append(distY / distance)
    return sin


def getErrors(cos, sin, residuals, radius):
    EaPlus = 0
    EaMinus = 0
    EbPlus = 0
    EbMinus = 0
    if len(cos) != len(sin):
        print("Error in getting (d1,d2): array of coordinates don't match!")
        return EaPlus, EbPlus, EaMinus, EbMinus
    if len(cos) != len(residuals):
        print("Error in getting (d1,d2): array of coordinates don't match!")
        return EaPlus, EbPlus, EaMinus, EbMinus
    for i in range(len(residuals)):
        diff = residuals[i] - radius
        if diff > 0:  # points outside circle
            EaPlus -= cos[i]
            EaMinus -= cos[i]
            EbPlus -= sin[i]
            EbMinus -= sin[i]
        if diff < 0:  # points inside circle
            EaPlus += cos[i]
            EaMinus += cos[i]
            EbPlus += sin[i]
            EbMinus += sin[i]
        if diff == 0:  # points on circle
            EaPlus += np.abs(cos[i])
            EaMinus -= np.abs(cos[i])
            EbPlus += np.abs(sin[i])
            EbMinus -= np.abs(sin[i])
    return EaPlus, EbPlus, EaMinus, EbMinus


def getAlphaBeta(EaPlus, EbPlus, EaMinus, EbMinus):
    alpha = 1 / 2
    beta = 1 / 2
    if EaPlus >= 0: alpha = 1
    if EaMinus <= 0: alpha = 0
    if EbPlus >= 0: beta = 1
    if EbMinus <= 0: beta = 0
    return alpha, beta
