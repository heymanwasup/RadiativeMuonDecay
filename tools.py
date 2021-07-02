import ROOT as R

#To-do
#   1.debug SetHistRange
#

class HistNameSvc:
    def __init__(self):
        self.Reset()

    def Reset(self):
        self.sampleTag = None
        self.channelTag = None
        self.selectTag = None

    def ResetPerEvent(self):
        self.channelTag = None
        self.selectTag = None

    def SetSampleTag(self,name):
        self.sampleTag = name

    def GetSampleTag(self):
        return self.sampleTag

    def SetChannelTag(self,name):
        self.channelTag = name
    
    def GetChannelTag(self):
        return self.channelTag

    def SetSelectTag(self,name):
        self.selectTag = name

    def GetSelectTag(self):
        return self.selectTag

    def GetFullName(self,vname):        
        tags = []
        if self.sampleTag!=None:
            tags.append(self.sampleTag)
        if self.channelTag!=None:
            tags.append(self.channelTag)
        if self.selectTag!=None:
            tags.append(self.selectTag)
        tags.append(vname)

        if len(tags)==1:
            full_name = tags[0]
        else:
            full_name = '_'.join(tags)
        return full_name

class HistSVc:
    def __init__(self):
        self.hists = {}
        self.name_handler = HistNameSvc()

    def BookFillTH1(self,vname,describ,n,x_i,x_f,value):
        name = self.name_handler.GetFullName(vname)
        if not name in self.hists:                    
            self.hists[name] = R.TH1F(name,describ,n,x_i,x_f)
        self.hists[name].Fill(value)

    def BookFillTH2(self,vname,describ,n_x,x_i,x_f,n_y,y_i,y_f,x,y):
        name = self.name_handler.GetFullName(vname)
        if not name in self.hists:
            self.hists[name] = R.TH2F(name,describ,n_x,x_i,x_f,n_y,y_i,y_f)
        self.hists[name].Fill(x,y)

    def BookFillTH3(self,vname,describ,n_x,x_i,x_f,n_y,y_i,y_f,n_z,z_i,z_f,x,y,z):
        name = self.name_handler.GetFullName(vname)
        if not name in self.hists:
            self.hists[name] = R.TH3F(name,describ,n_x,x_i,x_f,n_y,y_i,y_f,n_z,z_i,z_f)
        self.hists[name].Fill(x,y,z)

    def BookFullCutHist(self,name,describ,cuts_str,label):
        if not name in self.hists:
            self.hists[name] = R.TH1F(name,describ,len(cuts_str),-0.5,len(cuts_str)-0.5)
            for _nbin, _label in enumerate(cuts_str):
                self.hists[name].GetXaxis().SetBinLabel(_nbin+1,_label)
        self.hists[name].Fill(self.hists[name].GetXaxis().FindBin(label)-1)

    def GetHists(self):
        return self.hists

def Scale(hist):
    hist = hist.Clone()
    hist.Scale(1./hist.Integral())
    NbinsX  = hist.GetNbinsX()
    for nbin in range(NbinsX):
        hist.SetBinError(nbin+1,0.)
    return hist
        
def RetrieveHist(tfile,fmtStr):
    histA = tfile.Get(fmtStr.format('OMD')).Clone()
    histB = tfile.Get(fmtStr.format('RMD')).Clone()
    return histA,histB

def Draw(boxName,tfile,fmtStr,title,xTitle,yTitle,isScale=False,xRange=None):
    drawBox = DrawBox(boxName)
    drawBox.SetLeg('OMD','RMD',0.7,0.7,0.9,0.85)
    histA, histB = RetrieveHist(tfile,fmtStr)
    chi2ndf = histA.Chi2Test(histB,'WWCHI2/NDF')
    if isScale:
        histA, histB = Scale(histA),Scale(histB)
    drawBox.SetRange(xRange)
    canvas = drawBox.DrawIn(histA,histB,title,xTitle,yTitle,chi2ndf)
    return canvas,(histA,histB,drawBox)

class DrawBox:
    def __init__(self,name=''):
        self.xRange = None
        self.yRange = None
        self.name = name

    def SetLeg(self,nameA,nameB,x1,y1,x2,y2):
        self.legX1 = x1
        self.legY1 = y1
        self.legX2 = x2
        self.legY2 = y2
        self.nameA = nameA
        self.nameB = nameB

    def SetHist(self,hist,title,xTitle,yTitle,color):
        hist.SetLineColor(color)
        hist.SetStats(0)
        hist.SetTitle(title)
        hist.GetXaxis().SetTitle(xTitle)
        hist.GetYaxis().SetTitle(yTitle)
    def SetRange(self,xRange=None,yRange=None):
        self.xRange = xRange
        self.yRange = yRange

    def SetHistRange(self,histA,histB):
        if not self.yRange:
            maximum = max(histA.GetMaximum(),histB.GetMaximum())            
            histA.SetMaximum(1.1*maximum)
            histB.SetMaximum(1.1*maximum)
        else:
            minimum,maximum = tuple(self.yRange)
            histA.SetMaximum(1.1*maximum)
            histB.SetMaximum(1.1*maximum)
            histA.SetMinimum(minimum)
            histB.SetMinimum(minimum)

        if self.xRange!=None:
            histA.GetXaxis().SetRangeUser(*tuple(self.xRange))
            histB.GetXaxis().SetRangeUser(*tuple(self.xRange))

        
    def DrawIn(self,histA,histB,title,xTitle,yTitle,chi2ndf=-1.):
        self.SetHist(histA, title, xTitle, yTitle, 1)
        self.SetHist(histB, title, xTitle, yTitle, 2)
        self.SetHistRange(histA,histB)
        
        canvas = R.TCanvas(self.name+title,self.name+title)
        canvas.SetLeftMargin(0.2)
        self.leg = R.TLegend(self.legX1,self.legY1,self.legX2,self.legY2)
        self.leg.AddEntry(histA,self.nameA,'L')
        self.leg.AddEntry(histB,self.nameB,'L')
        canvas.cd(1)
        histA.Draw('same')
        histB.Draw('same')

        tex = R.TLatex()
#         tex.DrawLatex(0.7,0.64,'#color[4]{#Chi^{2}/ndf}: '+'{0:.1f}'.format(chi2ndf))
        text = '#Chi^{{2}}/ndf: {0:.1f}'.format(chi2ndf)
        tex.DrawLatexNDC(0.7,0.64,'#color[4]{' + text + '}')

        self.leg.Draw()

        return canvas

