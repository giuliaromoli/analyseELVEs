from src.utils.getPixelsInfos import getPixelStructure
import os
import numpy as np
import ROOT


def getPixelsInfosII(newPath, configFile, data, pixelsInfos, frameStart):
    if not os.path.exists(newPath):
        os.makedirs(newPath)

    pixelsInfos = []
    pixelsNRings = np.zeros((1, 20), dtype=int)
    filePixelInfo = open(newPath + "/infoPixels.txt", "w+")

    # get average of coherent pixels:
    allCandidates = [i for i in range(len(pixelsInfos)) if pixelsInfos[i][5] == 1]
    for time in range(frameStart, configFile.FrameEnd):
        candidates = [allCandidates[i] for i in range(len(allCandidates)) if (
                (pixelsInfos[allCandidates[i]][6] - 2 <= time) and (pixelsInfos[allCandidates[i]][6] + 2 >= time))]
        if len(candidates) > 0:
            plotPIXProfile = ROOT.TCanvas("plotPixelsProfile" + str(time))
            plotPIXP = ROOT.TMultiGraph()
            plotEneProfile = ROOT.TGraph()
            pixelProfile = []
            # plot pixels profile:
            for timeP in range(0, len(data)):
                sumEnergies = 0
                counter = 0
                for t in range(len(candidates)):
                    tt = candidates[t]
                    px = pixelsInfos[tt][1]
                    py = pixelsInfos[tt][2]
                    sumEnergies = sumEnergies + data[0, 0, timeP, px, py]
                    counter = counter + 1
                plotEneProfile.SetPoint(plotEneProfile.GetN(), timeP + frameStart, sumEnergies / counter)
                pixelProfile.append(sumEnergies / counter)
            plotEneProfile.SetMarkerStyle(20)
            plotEneProfile.SetMarkerColor(1)
            plotEneProfile.SetMarkerSize(0.5)
            plotPIXP.Add(plotEneProfile)
            plotPIXP.Draw("APL")
            plotPIXP.GetXaxis().SetTitle("Time")
            plotPIXP.GetYaxis().SetTitle("Energy [ADC]")
            plotPIXProfile.SaveAs(newPath + "/Profile" + str(time) + ".png")

    filePixelInfo.close()
