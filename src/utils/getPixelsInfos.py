import os
import ROOT
import numpy as np
import math
from scipy.optimize import curve_fit
import scipy.special as sse
from array import array
from src.utils.fillHisto import getPixelsPositions
from src.utils.fillHisto import fillHisto


def searchValley(data, limitThreshold):
    goldPrev = 0
    goldNext = 0
    maxEnergy = max(data)
    timeMax = np.where(data == maxEnergy)[0][0]
    # search for valley before:
    timeStart = timeMax - 1
    while timeStart > 0:
        goldPrev = timeStart - 1  # valley before
        if (goldPrev - 2) <= 0:
            break
        timePrev = timeStart - 2
        while (data[goldPrev] >= data[timePrev]) and (data[goldPrev] > 0):
            if (timePrev - 1) <= 0:
                break
            goldPrev = timePrev
            timePrev = timePrev - 1
        # rangeB=int(data[goldPrev]-nSigma*np.sqrt(data[goldPrev]))
        # rangeU=int(data[goldPrev]+nSigma*np.sqrt(data[goldPrev]))
        # ys=set(range(rangeB,rangeU))
        # yc=range(int(maxEnergy-nSigma*sigmaMaxEnergy),int(maxEnergy+nSigma*sigmaMaxEnergy))
        startRangeAboveValley = goldPrev - 5
        if startRangeAboveValley < 0:
            startRangeAboveValley = 0
        pointsAboveValley = [data[i] for i in range(goldPrev - 1, startRangeAboveValley, -1) if
                             data[i] >= data[goldPrev]]
        # print("valley before %d \n" %(goldPrev))
        maxPrev = goldPrev - 1  # local max before
        timePrev = goldPrev - 2
        while data[maxPrev] < data[timePrev]:
            if (timePrev - 1) <= 0:
                break
            maxPrev = timePrev
            timePrev = timePrev - 1
        # print("max before %d \n" %(maxPrev))
        if ((len(pointsAboveValley) < 4) and (
                data[goldPrev] - np.sqrt(data[goldPrev]) > limitThreshold)):  # (len(ys.intersection(yc))>0)and
            timeStart = timePrev
            continue
        if (data[maxPrev] - np.sqrt(data[maxPrev]) > limitThreshold) and (data[maxPrev] <= maxEnergy):
            if ((data[maxPrev - 1] + np.sqrt(data[maxPrev - 1])) <= (
                    data[goldPrev] - np.sqrt(data[goldPrev]))):
                timeStart = maxPrev
            else:
                break
        else:
            break
    # print(goldPrev)
    # search for next valley:
    timeStart = timeMax + 1
    while timeStart < len(data):
        goldNext = timeStart + 1  # valley after
        if (goldNext + 3) >= (len(data) - 1):
            break
        timeNext = timeStart + 2
        while (data[goldNext] >= data[timeNext]) and (data[goldNext] > 0):
            if (timeNext + 1) >= (len(data) - 2):
                break
            goldNext = timeNext
            timeNext = timeNext + 1
        # rangeB=int(data[goldNext]-np.sqrt(data[goldNext]))
        # rangeU=int(data[goldNext]+np.sqrt(data[goldNext]))
        # ys=set(range(rangeB,rangeU))
        # yc=range(int(maxEnergy-sigmaMaxEnergy),int(maxEnergy+sigmaMaxEnergy))
        endRangeAboveValley = goldNext + 5
        if endRangeAboveValley > len(data):
            endRangeAboveValley = len(data)
        pointsAboveValley = [data[i] for i in range(goldNext + 1, endRangeAboveValley) if
                             data[i] >= data[goldNext]]
        # print("valley next %d \n" %(goldNext))
        maxNext = goldNext + 1  # local max after
        timeNext = goldNext + 2
        while data[maxNext] < data[timeNext]:
            if (timeNext + 1) >= (len(data) - 2):
                break
            maxNext = timeNext
            timeNext = timeNext + 1
        # print("max next %d \n" %(maxNext))
        if ((len(pointsAboveValley) < 4) and (
                data[goldNext] - np.sqrt(data[goldNext]) >= limitThreshold)):  # (len(ys.intersection(yc))>0)and
            timeStart = timeNext
            continue
        if ((data[maxNext] - np.sqrt(data[maxNext]) > limitThreshold) and (
                data[maxNext] <= maxEnergy)):
            if ((data[maxNext + 1] + np.sqrt(data[maxNext + 1])) < (
                    data[goldNext] - np.sqrt(data[goldNext]))):
                timeStart = maxNext
            else:
                break
        else:
            break
    # print(goldNext)
    # print(goldPrev, goldNext)
    return goldPrev, goldNext


def gaussExpFit(x, y, x1, y1, sigma1, l1, threshold):
    l1 = 0.1
    area1 = y1 * sigma1 * 2
    try:
        (optParam, covMatrix) = curve_fit(
            lambda x, l1, sigma1, x1, area1: gaussExp(x, l1, sigma1, x1, area1, threshold), x, y,
            p0=[l1, sigma1, x1, area1])  # optParam optimal values, covMatrix covariance matrix
        # print(optParam)
        errParam = [covMatrix[0][0], covMatrix[1][1], covMatrix[2][2], covMatrix[3][3]]
        # print(errParam)
        if np.isinf(errParam[1]):
            optParam = [0, 0, 0, 0]
            errParam = [0, 0, 0, 0]
    except RuntimeError:
        optParam = [0, 0, 0, 0]
        errParam = [0, 0, 0, 0]
    return optParam, errParam


def gaussExp(x, l, s, m, area, threshold):
    # threshold=0
    gaussExp1 = 0.5 * l * np.exp(0.5 * l * (2 * m + l * np.power(s, 2) - 2 * x)) * sse.erfc(
        (m + l * np.power(s, 2) - x) / (np.sqrt(2) * s))  # exponential gaussian
    return area * gaussExp1 + threshold


def getPixelStructure(dataRaw, threshold, getSignal):

    data = dataRaw.copy()

    #for time in range(len(data)):
        #if time < 308:
            #data[time] = 0

    arrayPeaks = []
    optParamArray = []
    countZero = 0

    if getSignal == 1:
        canvasSignal = ROOT.TCanvas()
        signals = ROOT.TMultiGraph()

    while countZero < 5:  # is over threshold

        # print(countZero)

        if getSignal == 1:
            signal = ROOT.TGraph()
            for i in range(len(data)):
                signal.SetPoint(signal.GetN(), i, data[i])
            signal.SetMarkerStyle(20)
            signal.SetMarkerSize(0.5)
            signal.SetMarkerColor(countZero + 1)
            signals.Add(signal)

        maxEnergy = max(data)
        timeMax = np.where(data == maxEnergy)[0][0]
        sigmaMaxEnergy = np.sqrt(data[timeMax])

        if countZero == 0 and maxEnergy < (threshold + np.sqrt(threshold)):
            break

        if countZero > 0 and maxEnergy < 5:
            break

        if countZero == 0 and maxEnergy < 30:
            break

        thresholdSinglePeak = threshold
        if countZero > 0:
            thresholdSinglePeak = 0

        (goldPrev, goldNext) = searchValley(data, thresholdSinglePeak)

        # set initial conditions for fitting:
        sigma1 = int((timeMax - goldPrev) / 2)
        lambda1 = np.log(2) / (np.sqrt(2 * np.log(2)) * sigma1)
        intTime1 = goldPrev
        if intTime1 < 0:
            intTime1 = 0
        intTime2 = goldNext
        if intTime2 > len(data) - 1:
            intTime2 = len(data) - 1
        # print(range(intTime1,intTime2))

        if countZero == 0:
            dataToAverage = [data[i] for i in range(0, intTime1) if data[i] > 0]
            thresholdSinglePeak = np.average(dataToAverage)
            threshold = thresholdSinglePeak

        if len(data[intTime1:intTime2]) > 4:  # at least 5 points to fit
            (optParam, errParam) = gaussExpFit(np.array(range(intTime1, intTime2)), data[intTime1:intTime2], timeMax,
                                               maxEnergy, sigma1, lambda1, thresholdSinglePeak)
            # print(optParam)
            if (optParam[2] > 0) and (optParam[3] > 0) and (optParam[3] < 1000):
                fitCurve = gaussExp(np.array(range(len(data))), optParam[0], optParam[1], optParam[2], optParam[3],
                                    thresholdSinglePeak)
                fitCurvePos = [fitCurve[i] for i in range(len(fitCurve)) if fitCurve[i] >= 0.1]
                if len(fitCurvePos) > 0:
                    maxEnergyFit = max(fitCurvePos)
                    if maxEnergyFit > 0:
                        # get infos:
                        fullWidth = goldNext - goldPrev
                        deExcTime = goldNext - timeMax
                        totEnergy = np.sum(data[goldPrev:goldNext])
                        sigmaTotEnergy = np.average(np.sqrt(data[goldPrev:goldNext]))
                        halfPeak = maxEnergy / 2
                        abovePeak = [data[i] for i in range(len(data)) if data[i] >= halfPeak]
                        energyPeak = np.average(abovePeak)
                        sigmaEnergyPeak = np.average(np.sqrt(abovePeak))
                        sigmaMaxEnergyFit = np.sqrt(maxEnergyFit)
                        halfPeakFit = maxEnergyFit / 2
                        abovePeakFit = [fitCurve[i] for i in range(len(fitCurve)) if fitCurve[i] >= halfPeakFit]
                        energyPeakFit = np.average(abovePeakFit)
                        sigmaEnergyPeakFit = np.average(np.sqrt(abovePeakFit))
                        arrayPeaks.append((timeMax, maxEnergy, sigmaMaxEnergy, fullWidth, deExcTime, totEnergy,
                                           sigmaTotEnergy, energyPeak, sigmaEnergyPeak, goldPrev, goldNext,
                                           thresholdSinglePeak))
                        optParamArray.append(
                            (optParam, maxEnergyFit, sigmaMaxEnergyFit, energyPeakFit, sigmaEnergyPeakFit))
                        # subtract fit to initial signal:
                        dataNew = data
                        for time in range(len(data)):
                            if goldPrev <= time <= goldNext:
                                dataNew[time] = 0
                                continue
                            yToSub = gaussExp(time, optParam[0], optParam[1], optParam[2], optParam[3],
                                              thresholdSinglePeak)
                            if (yToSub < 0.1) or (math.isnan(yToSub)):
                                yToSub = thresholdSinglePeak
                            # print(time,yToSub)
                            if yToSub > 0:
                                dataNew[time] = data[time] - yToSub
                                if dataNew[time] <= 0:
                                    dataNew[time] = 0
                            else:
                                dataNew[time] = data[time]
                            # else:
                            # dataNew[time]=0
                        data = dataNew
                    else:
                        break
                else:
                    break
            else:
                fullWidth = goldNext - goldPrev
                deExcTime = goldNext - timeMax
                totEnergy = np.sum(data[goldPrev:goldNext])
                sigmaTotEnergy = np.average(np.sqrt(data[goldPrev:goldNext]))
                halfPeak = maxEnergy / 2
                abovePeak = [data[i] for i in range(len(data)) if data[i] >= halfPeak]
                energyPeak = np.average(abovePeak)
                sigmaEnergyPeak = np.average(np.sqrt(abovePeak))
                arrayPeaks.append((timeMax, maxEnergy, sigmaMaxEnergy, fullWidth, deExcTime, totEnergy,
                                   sigmaTotEnergy, energyPeak, sigmaEnergyPeak, goldPrev, goldNext,
                                   thresholdSinglePeak))
                optParamArray.append(([0, 0, 0, 0], 0, 0, 0, 0))
                break
        else:
            break
        countZero = countZero + 1

    # print(arrayPeaks)
    # print(optParamArray)

    if getSignal == 1:
        signals.Draw("APL")
        canvasSignal.SaveAs(r'/home/jule/Scrivania/Mini-Euso/fitCircle/signal.png')

    return arrayPeaks, optParamArray


def getPixelsInfos(newPath, configFile, data, frameStart):
    if not os.path.exists(newPath):
        os.makedirs(newPath)
    if not os.path.exists(newPath + '/badPixels'):
        os.makedirs(newPath + '/badPixels')
    pixelsInfos = []
    pixelsNRings = np.zeros((1, 20), dtype=int)
    filePixelInfo = open(newPath + "/infoPixels.txt", "w+")
    # analyse pixel:
    (kmdistanceX, kmdistanceY) = getPixelsPositions(configFile.PixelFileName, configFile.HH)
    for px in range(48):
        for py in range(48):
            goodPlot = 0
            pixelLabel = px + 48 * py  # pixel label
            pixX = px  # x-label
            pixY = py  # y-label
            xCoord = -kmdistanceX[px][py]  # x-coord
            yCoord = -kmdistanceY[px][py]  # y-coord
            nPeak = 0
            maxEnergy = max(data[0, 0, :, px, py])
            if ((configFile.limXmin <= xCoord <= configFile.limXmax) and (
                    configFile.limYmin <= yCoord <= configFile.limYmax) and (
                    maxEnergy - np.sqrt(maxEnergy) > configFile.Threshold + np.sqrt(configFile.Threshold))):
                # print pixel profile:
                plotPIXProfile = ROOT.TCanvas("plotPixProfile" + str(pixelLabel))
                plotPIX = ROOT.TMultiGraph()
                arrayX = np.zeros(len(data[0, 0, :, px, py]), dtype=int)
                arrayY = np.zeros(len(data[0, 0, :, px, py]), dtype=float)
                arrayErrY = np.zeros(len(data[0, 0, :, px, py]), dtype=float)
                zeros = np.zeros(len(data[0, 0, :, px, py]), dtype=int)
                for time in range(len(data[0, 0, :, px, py])):
                    arrayX[time] = time + frameStart
                    arrayY[time] = data[0, 0, time, px, py]
                    arrayErrY[time] = np.sqrt(data[0, 0, time, px, py])
                plotEnePixelErr = ROOT.TGraphErrors(len(arrayX), array('f', arrayX), array('f', arrayY),
                                                    array('f', zeros), array('f', arrayErrY))
                plotEnePixelErr.SetMarkerStyle(20)
                plotEnePixelErr.SetMarkerColor(1)
                plotEnePixelErr.SetLineColor(27)
                plotEnePixelErr.SetMarkerSize(0.5)
                plotPIX.Add(plotEnePixelErr)
                plotEnePixel = ROOT.TGraph()
                for time in range(len(data[0, 0, :, px, py])):
                    plotEnePixel.SetPoint(plotEnePixel.GetN(), time + frameStart,
                                          data[0, 0, time, px, py])
                plotEnePixel.SetMarkerStyle(20)
                plotEnePixel.SetMarkerColor(1)
                plotEnePixel.SetMarkerSize(0.5)
                plotPIX.Add(plotEnePixel)
                # get infos:
                (arrayPeaks, optParamArray) = getPixelStructure(
                    data[0, 0, 0:(configFile.FrameEnd - frameStart), px, py], configFile.Threshold, 0)
                if len(arrayPeaks) > 0:
                    thresholdSinglePeak = arrayPeaks[0][11]
                    # get timeStartPeak:
                    if len(arrayPeaks) > 1:
                        secondP = 1
                        while abs(arrayPeaks[0][0] - arrayPeaks[secondP][0]) > 20:
                            secondP = secondP + 1
                            if secondP == len(arrayPeaks):
                                if abs(arrayPeaks[0][0] - arrayPeaks[secondP - 1][0]) > 20:
                                    secondP = 0
                                break
                        # print(secondP)
                        if arrayPeaks[0][0] < arrayPeaks[secondP][0]:
                            timeStartPeak = arrayPeaks[0][0] - (arrayPeaks[0][3] - arrayPeaks[0][4]) + frameStart
                        else:
                            timeStartPeak = arrayPeaks[secondP][0] - (
                                    arrayPeaks[secondP][3] - arrayPeaks[secondP][4]) + frameStart
                    else:
                        timeStartPeak = arrayPeaks[0][0] - (arrayPeaks[0][3] - arrayPeaks[0][4]) + frameStart
                    timeStart = timeStartPeak - frameStart
                    while timeStart > 0 and data[0, 0, timeStart, px, py] > thresholdSinglePeak:
                        timeStart = timeStart - 1
                    if timeStart + frameStart < timeStartPeak:  # goldPrev is changed
                        arrayPeaksList = []
                        for s in range(len(arrayPeaks)):
                            arrayPeaksList.append(list(arrayPeaks[s][:]))
                        if len(arrayPeaks) > 1:
                            if arrayPeaks[0][0] < arrayPeaks[secondP][0]:
                                arrayPeaksList[0][9] = timeStart  # goldPrev
                                arrayPeaksList[0][3] = arrayPeaksList[0][10] - timeStart  # fullWidth
                                arrayPeaksList[0][5] = np.sum(
                                    data[0, 0, arrayPeaksList[0][9]:arrayPeaksList[0][10], px, py])  # totEnergy
                                arrayPeaksList[0][6] = np.average(np.sqrt(
                                    data[0, 0, arrayPeaksList[0][9]:arrayPeaksList[0][10], px, py]))  # sigmaTotEnergy
                            else:
                                arrayPeaksList[secondP][9] = timeStart  # goldPrev
                                arrayPeaksList[secondP][3] = arrayPeaksList[secondP][10] - timeStart  # fullWidth
                                arrayPeaksList[secondP][5] = np.sum(
                                    data[0, 0, arrayPeaksList[secondP][9]:arrayPeaksList[secondP][10], px,
                                    py])  # totEnergy
                                arrayPeaksList[secondP][6] = np.average(np.sqrt(
                                    data[0, 0, arrayPeaksList[secondP][9]:arrayPeaksList[secondP][10], px,
                                    py]))  # sigmaTotEnergy
                        else:
                            arrayPeaksList[0][9] = timeStart  # goldPrev
                            arrayPeaksList[0][3] = arrayPeaksList[0][10] - timeStart  # fullWidth
                            arrayPeaksList[0][5] = np.sum(
                                data[0, 0, arrayPeaksList[0][9]:arrayPeaksList[0][10], px, py])  # totEnergy
                            arrayPeaksList[0][6] = np.average(np.sqrt(
                                data[0, 0, arrayPeaksList[0][9]:arrayPeaksList[0][10], px, py]))  # sigmaTotEnergy
                        arrayPeaks = []
                        for s in range(len(arrayPeaksList)):
                            arrayPeaks.append(tuple(arrayPeaksList[s][:]))
                    timeStartPeak = timeStart + frameStart
                    # get timeEndPeak:
                    if len(arrayPeaks) > 1:
                        if arrayPeaks[0][0] < arrayPeaks[secondP][0]:
                            timeEndPeak = arrayPeaks[secondP][0] + arrayPeaks[secondP][4] + frameStart
                        else:
                            timeEndPeak = arrayPeaks[0][0] + arrayPeaks[0][4] + frameStart
                    else:
                        timeEndPeak = arrayPeaks[0][0] + arrayPeaks[0][4] + frameStart
                    timeEnd = timeEndPeak - frameStart
                    while timeEnd < len(data) and data[0, 0, timeEnd, px, py] > thresholdSinglePeak:
                        timeEnd = timeEnd + 1
                    if timeEnd + frameStart > timeEndPeak:  # goldNext is changed
                        arrayPeaksList = []
                        for s in range(len(arrayPeaks)):
                            arrayPeaksList.append(list(arrayPeaks[s][:]))
                        if len(arrayPeaks) > 1:
                            if arrayPeaks[0][0] < arrayPeaks[secondP][0]:
                                arrayPeaksList[0][10] = timeEnd  # goldNext
                                arrayPeaksList[0][4] = timeEnd - arrayPeaksList[0][0]  # deExcTime
                                arrayPeaksList[0][3] = timeEnd - arrayPeaksList[0][9]  # fullWidth
                                arrayPeaksList[0][5] = np.sum(
                                    data[0, 0, arrayPeaksList[0][9]:arrayPeaksList[0][10], px, py])  # totEnergy
                                arrayPeaksList[0][6] = np.average(np.sqrt(
                                    data[0, 0, arrayPeaksList[0][9]:arrayPeaksList[0][10], px, py]))  # sigmaTotEnergy
                            else:
                                arrayPeaksList[secondP][10] = timeEnd  # goldNext
                                arrayPeaksList[secondP][4] = timeEnd - arrayPeaksList[secondP][0]  # deExcTime
                                arrayPeaksList[secondP][3] = timeEnd - arrayPeaksList[secondP][9]  # fullWidth
                                arrayPeaksList[secondP][5] = np.sum(
                                    data[0, 0, arrayPeaksList[secondP][9]:arrayPeaksList[secondP][10], px,
                                    py])  # totEnergy
                                arrayPeaksList[secondP][6] = np.average(np.sqrt(
                                    data[0, 0, arrayPeaksList[secondP][9]:arrayPeaksList[secondP][10], px,
                                    py]))  # sigmaTotEnergy
                        else:
                            arrayPeaksList[0][10] = timeEnd  # goldNext
                            arrayPeaksList[0][4] = timeEnd - arrayPeaksList[0][0]  # deExcTime
                            arrayPeaksList[0][3] = timeEnd - arrayPeaksList[0][9]  # fullWidth
                            arrayPeaksList[0][5] = np.sum(
                                data[0, 0, arrayPeaksList[0][9]:arrayPeaksList[0][10], px, py])  # totEnergy
                            arrayPeaksList[0][6] = np.average(np.sqrt(
                                data[0, 0, arrayPeaksList[0][9]:arrayPeaksList[0][10], px, py]))  # sigmaTotEnergy
                        arrayPeaks = []
                        for s in range(len(arrayPeaksList)):
                            arrayPeaks.append(tuple(arrayPeaksList[s][:]))
                    timeEndPeak = timeEnd + frameStart
                    totDuration = timeEndPeak - timeStartPeak
                    totDeExcitation = timeEndPeak - int(arrayPeaks[0][0]) - frameStart
                    for p in range(len(arrayPeaks)):
                        nPeak = nPeak + 1
                        pixelsNRings[0][nPeak] = pixelsNRings[0][nPeak] + 1
                        goodPlot = 1
                        timeMax = int(arrayPeaks[p][0]) + frameStart
                        maxPEnergy = arrayPeaks[p][1]
                        sigmaMaxEnergy = arrayPeaks[p][2]
                        fullWidth = int(arrayPeaks[p][3])
                        deExcTime = int(arrayPeaks[p][4])
                        totEnergy = arrayPeaks[p][5]
                        sigmaTotEnergy = arrayPeaks[p][6]
                        thresholdSinglePeak = arrayPeaks[p][11]
                        energyPeak = arrayPeaks[p][7]
                        sigmaEnergyPeak = arrayPeaks[p][8]
                        if optParamArray[p][0][0] > 0:
                            tauFit = 1 / optParamArray[p][0][0]
                            sigmaFit = np.sqrt(np.power(optParamArray[p][0][1], 2) + np.power(tauFit, 2))
                            meanFit = int(optParamArray[p][0][2]) + tauFit + frameStart
                            areaFit = optParamArray[p][0][3]
                            maxEnergyFit = optParamArray[p][1]
                            sigmaMaxEnergyFit = optParamArray[p][2]
                            energyPeakFit = optParamArray[p][3]
                            sigmaEnergyPeakFit = optParamArray[p][4]
                            kFit = np.power(tauFit / sigmaFit, 3)
                            lambdaFit = np.power(kFit, 1 / 3)
                        else:
                            tauFit = 0
                            sigmaFit = 0
                            meanFit = 0
                            areaFit = 0
                            maxEnergyFit = 0
                            sigmaMaxEnergyFit = 0
                            energyPeakFit = 0
                            sigmaEnergyPeakFit = 0
                            kFit = 0
                            lambdaFit = 0
                        filePixelInfo.write(
                            "%d \t %d \t %d \t %5.2f \t %5.2f \t %d \t %d \t %5.2f \t %5.2f \t %d \t %d \t %5.2f \t "
                            "%5.2f \t %5.2f \t %5.2f \t %5.2f \t %5.2f \t %d \t %5.2f \t %5.2f \t %5.2f \t %5.2f \t "
                            "%5.2f \t %5.2f \t %5.2f \t %d \t %d \n" % (
                                pixelLabel, pixX, pixY, xCoord, yCoord, nPeak, timeMax, maxPEnergy, sigmaMaxEnergy,
                                fullWidth, deExcTime, totEnergy, sigmaTotEnergy, energyPeak, sigmaEnergyPeak, tauFit,
                                sigmaFit, meanFit, areaFit, maxEnergyFit, sigmaMaxEnergyFit, energyPeakFit,
                                sigmaEnergyPeakFit, thresholdSinglePeak, lambdaFit, totDuration, totDeExcitation))
                fitTot = ROOT.TGraph()
                xToPlot = np.array(range(frameStart, frameStart + len(data[0, 0, :, px, py])))
                fitTotArray = np.zeros((1, len(xToPlot)))
                for oo in range(len(optParamArray)):
                    if optParamArray[oo][0][0] > 0:
                        yToPlot = gaussExp(xToPlot, optParamArray[oo][0][0], optParamArray[oo][0][1],
                                           optParamArray[oo][0][2] + frameStart, optParamArray[oo][0][3],
                                           arrayPeaks[oo][11])
                        if oo > 0:
                            yToPlot = gaussExp(xToPlot, optParamArray[oo][0][0], optParamArray[oo][0][1],
                                               optParamArray[oo][0][2] + frameStart, optParamArray[oo][0][3],
                                               arrayPeaks[oo][11])
                        fitExp = ROOT.TGraph()
                        for xx in range(len(xToPlot)):
                            if yToPlot[xx] >= 0:
                                fitTotArray[0][xx] = fitTotArray[0][xx] + yToPlot[xx]
                                fitExp.SetPoint(fitExp.GetN(), xToPlot[xx], yToPlot[xx])
                                fitExp.SetLineStyle(2)
                        fitExp.SetLineColor(oo + 91)
                        fitExp.SetLineWidth(3)
                        plotPIX.Add(fitExp)
                for xx in range(len(xToPlot)):
                    fitTot.SetPoint(fitTot.GetN(), xToPlot[xx], fitTotArray[0][xx])
                fitTot.SetLineStyle(2)
                fitTot.SetLineColor(4)
                fitTot.SetLineWidth(3)
                plotPIX.Add(fitTot)
                plotPIX.Draw("APL")
                if len(arrayPeaks) > 0:
                    line1 = ROOT.TLine(timeStartPeak, 0, timeStartPeak, 1000)
                    line1.SetLineColor(2)
                    line1.SetLineStyle(9)
                    line1.DrawLine(timeStartPeak, 0, timeStartPeak, 1000)
                    line2 = ROOT.TLine(timeEndPeak, 0, timeEndPeak, 1000)
                    line2.SetLineColor(2)
                    line2.SetLineStyle(9)
                    line2.DrawLine(timeEndPeak, 0, timeEndPeak, 1000)
                    line3 = ROOT.TLine(timeEndPeak, 0, timeEndPeak, 1000)
                    line3.SetLineColor(4)
                    line3.SetLineStyle(9)
                    line3.DrawLine(frameStart, arrayPeaks[0][11], frameStart + 450, arrayPeaks[0][11])

                plotPIX.GetXaxis().SetTitle("Time")
                plotPIX.GetXaxis().SetRangeUser(frameStart, configFile.FrameEnd)
                plotPIX.GetYaxis().SetTitle("Energy [ADC]")
                plotPIX.GetXaxis().SetRangeUser(960, 1040)
                if goodPlot == 1:
                    plotPIXProfile.SaveAs(newPath + "/" + str(pixelLabel) + ".png")
                else:
                    pixelsNRings[0][0] = pixelsNRings[0][0] + 1
                    plotPIXProfile.SaveAs(newPath + "/badPixels/" + str(pixelLabel) + ".png")
    filePixelInfo.close()
    print(pixelsNRings)
    return pixelsInfos


def getPixelsPlots(newPath, pixelsInfos, configFile):
    if not os.path.exists(newPath):
        os.makedirs(newPath)
    filePixelInfo = open(newPath + "/infoPixels.txt", "w+")
    # plots of main peak:
    candidates = [i for i in range(len(pixelsInfos)) if (pixelsInfos[i][5] == 1)]

    # secondPeakCandidates = []
    # for t in range(len(candidates)):
    #     tt = candidates[t]
    #     timeMax = pixelsInfos[tt][6]
    #     peaks = [x for x in range(len(pixelsInfos)) if (pixelsInfos[x][0] == pixelsInfos[tt][0]) and (pixelsInfos[x][5] > 1) and (pixelsInfos[x][6]-timeMax > 0)]
    #     if len(peaks) > 0:
    #         diffs = []
    #         for p in range(len(peaks)):
    #             diffs.append(abs(timeMax-pixelsInfos[peaks[p]][6]))
    #         minDiffs=min(diffs)
    #         secondP = [x for x in range(len(peaks)) if diffs[x]== minDiffs][0]
    #         secondPeakCandidates.append(peaks[secondP])
    # candidates = secondPeakCandidates

    # tot-duration:
    plotPixelsInfos(newPath, candidates, pixelsInfos, configFile, "totDuration", 25, 0, 50, filePixelInfo)
    # tot-de-excitation:
    plotPixelsInfos(newPath, candidates, pixelsInfos, configFile, "totDe-Excitation", 26, 0, 30, filePixelInfo)
    # duration:
    #plotPixelsInfos(newPath, candidates, pixelsInfos, configFile, "duration", 9, 0, 50, filePixelInfo)
    # de-excitation:
    #plotPixelsInfos(newPath, candidates, pixelsInfos, configFile, "de-Excitation", 10, 0, 30, filePixelInfo)
    # threshold:
    plotPixelsInfos(newPath, candidates, pixelsInfos, configFile, "background", 23, 0, 15, filePixelInfo)
    # sigma:
    plotPixelsInfos(newPath, candidates, pixelsInfos, configFile, "sigma", 16, 0, 30, filePixelInfo)
    # tau:
    plotPixelsInfos(newPath, candidates, pixelsInfos, configFile, "tau", 15, 0, 30, filePixelInfo)
    # lambda:
    plotPixelsInfos(newPath, candidates, pixelsInfos, configFile, "lambda", 24, 0, 1, filePixelInfo)
    # maxEnergy:
    plotPixelsInfos(newPath, candidates, pixelsInfos, configFile, "maxEnergy", 7, 0, 120, filePixelInfo)
    # time distance:
    # deltaT:
    # pointsToPlot = np.zeros((48, 48), dtype=np.float32)
    # ave = []
    # for t in range(len(candidates)):
    #     tt = candidates[t]
    #     timeMax = pixelsInfos[tt][6]
    #     peaks = [x for x in range(len(pixelsInfos)) if
    #              (pixelsInfos[x][0] == pixelsInfos[tt][0]) and (pixelsInfos[x][5] > 1) and (pixelsInfos[x][6]-timeMax > 0)]
    #     if len(peaks) > 0:
    #         diffs = []
    #         for p in range(len(peaks)):
    #             diffs.append(abs(timeMax - pixelsInfos[peaks[p]][6]))
    #         minDiffs = min(diffs)
    #         pixX = pixelsInfos[tt][1]
    #         pixY = pixelsInfos[tt][2]
    #         pointsToPlot[pixX][pixY] = minDiffs  # delta T
    #         ave.append(minDiffs)
    # print("deltaT: %5.2f " % (np.average(ave)))
    # print("error: %5.2f " % (np.std(ave)))
    # canvasH = ROOT.TCanvas()
    # histoH = ROOT.TH2Poly()
    # (histoH, canvasH) = fillHisto(newPath, configFile.PixelFileName, configFile.HH, pointsToPlot, 0, 35, [], [], 1, "", "/deltaT.png")
    # histoH.Draw("COLZ")
    # canvasH.Update()
    # canvasH.SaveAs(newPath + "/deltaT.png")
    #
    # canvasH = ROOT.TCanvas()
    # histoDuration = ROOT.TH1F("TH", "", 35, 0, 35)
    # for t in range(len(ave)):
    #     histoDuration.Fill(ave[t], 1)
    # histoDuration.SetStats(0)
    # histoDuration.SetFillColor(38)
    # histoDuration.SetFillStyle(3001)
    # histoDuration.Draw()
    # histoDuration.GetXaxis().SetTitle("time distance [GTU]")
    # histoDuration.GetYaxis().SetTitle("counts")
    # canvasH.SaveAs(newPath + "/histoDeltaT.png")

    filePixelInfo.close()


def plotPixelsInfos(newPath, candidates, pixelsInfos, configFile, propName, propLabel, yMin, yMax, filePixelInfo):
    pointsToPlot = np.zeros((48, 48), dtype=np.float32)
    ave = []
    for t in range(len(candidates)):
        tt = candidates[t]
        pixX = pixelsInfos[tt][1]
        pixY = pixelsInfos[tt][2]
        pointsToPlot[pixX][pixY] = pixelsInfos[tt][propLabel]
        ave.append(pixelsInfos[tt][propLabel])
        # if (pixelsInfos[tt][26]>0):
        #     print(pixelsInfos[tt][0])
    filePixelInfo.write(propName + ": %5.2f " % (np.average(ave)) + "\n")
    filePixelInfo.write("error: %5.2f " % (np.std(ave)) + "\n")
    canvasH = ROOT.TCanvas()
    histoH = ROOT.TH2Poly()
    (histoH, canvasH) = fillHisto(newPath, configFile.PixelFileName, configFile.HH, pointsToPlot, yMin, yMax, [], [], 1,
                                  "",
                                  "/" + propName + ".png")
    histoH.Draw("COLZ")
    canvasH.Update()
    canvasH.SaveAs(newPath + "/" + propName + ".png")

    canvasH = ROOT.TCanvas()
    histoDuration = ROOT.TH1F("TD", "", 25, yMin, yMax)
    for t in range(len(ave)):
        histoDuration.Fill(ave[t], 1)
    histoDuration.SetStats(0)
    histoDuration.SetFillColor(38)
    histoDuration.SetFillStyle(3001)
    histoDuration.Draw()
    histoDuration.GetXaxis().SetTitle(propName)
    histoDuration.GetYaxis().SetTitle("counts")
    canvasH.SaveAs(newPath + "/histo" + propName + ".png")


def getPhiEnergy(newPath, pixelsInfos, ElveCentreX, ElveCentreY):
    # this function returns phi angles for selected points:
    # (time,posX,posY,energy)->(time,phi,energy)

    pixelsInfosPhi = []
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


if __name__ == '__main__':
    pass
