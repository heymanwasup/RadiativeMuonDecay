#read OMD/RMD trees, write histos into root file

import ROOT as R
import tools

R.gROOT.SetBatch(True)

class Event:
    def __init__(self):
        self.reset()

    def reset(self):
        self.Etotal = 0.
        
        self.E1 = -1.
        self.E2 = -1.
        self.E3 = -1.


        self.Et1 = -1.
        self.Et2 = -1.
        self.Et3 = -1.

        self.N1 = -1
        self.N2 = -1
        self.N3 = -1

        self.Nt1 = -1
        self.Nt2 = -1
        self.Nt3 = -1

        #ns
        self.T1 = -999.
        self.T2 = -999.
        self.T3 = -999.

        self.Tt1 = -999.
        self.Tt2 = -999.
        self.Tt3 = -999.

        self.caloHitX = -999.
        self.caloHitY = -999.
        self.caloHitZ = -999.

        self.nCount200 = 0 #num of calos with E>200
        self.nCount100 = 0 #num of calos with E>100
        self.nCount50 = 0 #num of calos with E>100
        self.nCount0   = 0 #num of calos with E>0
        self.posiTrackID = -1
        self.time_energy_dict = {}
        self.time_ncalo_dict = {}

class ReadVariables:
    def __init__(self):
        # self.tfile_omd = R.TFile('./input/gm2offline_ana_omd.root')
        # self.tfile_rmd = R.TFile('./input/gm2offline_ana_rmd.root')
        self.tfile_omd = R.TFile('./input/gm2offline_ana_omd_all.root')
        self.tfile_rmd = R.TFile('./input/gm2offline_ana_rmd_all.root')
        self.ttree_omd = self.tfile_omd.Get('G2PhaseAnalyzer/g2phase')
        self.ttree_rmd = self.tfile_rmd.Get('G2PhaseAnalyzer/g2phase')

        self.hist_svc = tools.HistSVc()        
        self.name_handler = self.hist_svc.name_handler

        self.counter = [
            'Total',

            'OMD',
            'OMD_omd',
            'OMD_rmd',            

            'RMD',
            'RMD_omd',
            'RMD_rmd',            
        ]
        

    def WriteHistos(self):
        tfile_out = R.TFile('./output/hists_output_all_v1.root','recreate')
        tfile_out.cd()

        self.loopTree(self.ttree_omd,'OMD')
        self.loopTree(self.ttree_rmd,'RMD')

        tfile_out.cd()
        tfile_out.Write()
        tfile_out.Close()

    def fillEvent(self,ttree,event):
        
        event.posiTrackID = ttree.posiTrackID
        event.isCaloHit = ttree.isCaloHit
        if event.isCaloHit:
            event.caloHitX = ttree.caloHitX
            event.caloHitY = ttree.caloHitY
            event.caloHitZ = ttree.caloHitZ

        time_energy_dict = event.time_energy_dict
        time_ncalo_dict = event.time_ncalo_dict
        for nCalo in range(24):
            E = ttree.caloHitEnergyArray[nCalo]
            T = ttree.caloHitTimeArray[nCalo]
            if E==0. and T==0:
                continue
            time_energy_dict[T] = E
            time_ncalo_dict[T] = nCalo

            event.Etotal += E

            if E > 200.:
                event.nCount200 += 1
            if E > 100.:
                event.nCount100 += 1
            if E > 50.:
                event.nCount50 += 1
            if E > 0.:
                event.nCount0 += 1

            if E > event.E1:
                event.E1 = E
                event.N1 = nCalo
                event.T1 = T
            elif E > event.E2:
                event.E2 = E
                event.N2 = nCalo
                event.T2 = T
            elif E > event.E3:
                event.E3 = E
                event.N3 = nCalo
                event.T3 = T

        NCalo = len(time_energy_dict)
        time_list = time_energy_dict.keys()
        time_list.sort()

        if NCalo>0:
            event.Et1 = time_energy_dict[time_list[0]]
            event.Tt1 = time_list[0]
            event.Nt1 = time_ncalo_dict[time_list[0]]
        if NCalo>1:
            event.Et2 = time_energy_dict[time_list[1]]
            event.Tt2 = time_list[1]
            event.Nt2 = time_ncalo_dict[time_list[1]]

        if NCalo>2:
            event.Et3 = time_energy_dict[time_list[2]]
            event.Tt3 = time_list[2]
            event.Nt3 = time_ncalo_dict[time_list[2]]
        
        return True


    def loopTree(self,ttree,treeType):

        self.name_handler.Reset()
        self.name_handler.SetSampleTag(treeType)
        
        event = Event()
        Nentries = ttree.GetEntries()    
        for n in range(Nentries):
            if not n%50000:
                print 'Processed %s events for %s tree'%(n,treeType)
            ttree.GetEntry(n)            
           
            event.reset()
            self.name_handler.ResetPerEvent()
            
            self.fillEvent(ttree,event)
            self.fillHists(event)

    def fillAnalysisHists(self,event):
                
        self.hist_svc.BookFillTH1(
            'Etotal','total calo energy deposit',50,0,3000, event.Etotal)

        self.hist_svc.BookFillTH1(
            'nCaloE200', 'num of calos with E>200 MeV', 5, -0.5, 4.5, event.nCount200)
        self.hist_svc.BookFillTH1(
            'nCaloE100', 'num of calos with E>100 MeV', 5, -0.5, 4.5, event.nCount100)
        self.hist_svc.BookFillTH1(
            'nCaloE50', 'num of calos with E>50 MeV', 5, -0.5, 4.5, event.nCount50)
        self.hist_svc.BookFillTH1(
            'nCaloE0', 'num of calos with E>0 MeV', 5, -0.5, 4.5, event.nCount0)



        if event.Et1 > 0:
            self.hist_svc.BookFillTH1(
                'Et1','1st coming Xtal Edep',50,0,3000,event.Et1)
        if event.Et2 > 0:
            self.hist_svc.BookFillTH1(
                'Et2','2nd coming Xtal Edep',50,0,3000,event.Et2)
        if event.Et3 > 0:
            self.hist_svc.BookFillTH1(
                'Et3','3rd coming Xtal Edep',50,0,3000,event.Et3)

        if event.E1 > 0:
            self.hist_svc.BookFillTH1(
                'E1','1st calo energy deposit',50,0,3000,event.E1)
        if event.E2 > 0:
            self.hist_svc.BookFillTH1(
                'E2','2nd calo energy deposit',50,0,500,event.E2)
        if event.E3 > 0:
            self.hist_svc.BookFillTH1(
                'E3','3rd calo energy deposit',50,0,200,event.E3)

        self.hist_svc.BookFillTH2(
            'Et1E1','Xtal E dep.: earliest vs highest',50,0,3000,50,0,3000,event.Et1,event.E1)

        if event.E2>0:
            self.hist_svc.BookFillTH2(
                'E1E2','Xtal E dep.: leading vs sub-leading',50,0,3000,50,0,3000,event.E1,event.E2)

            self.hist_svc.BookFillTH2(
                'Et1Et2','Xtal E dep.: first vs second',50,0,3000,50,0,3000,event.Et1,event.Et2)

        deltaN = min(24-abs(event.N2-event.N1),abs(event.N2-event.N1))
        deltaT = abs(event.T1 - event.T2)

        # deltaNt = min(24-abs(event.Nt2-event.Nt1),abs(event.Nt2-event.Nt1))
        # deltaTt = abs(event.Tt1 - event.Tt2)

        count_E = zip([event.nCount200,event.nCount100,event.nCount50,event.nCount0],[200,100,50,0])
        for nCount,E in count_E:
            #deltaN for leading and sub-leading calos with E>200,100,50,0
            if nCount >= 2:                        
                self.hist_svc.BookFillTH1('E%d_deltaN'%(E),'Delta N of calos with E>%d'%(E),14,-0.5,13.5,deltaN)
                #deltaT for nearby (in space) 2 calos with E>200,100,50,0
                if deltaN == 1:
                    self.hist_svc.BookFillTH1('E%d_deltaT'%(E),'dT of nearby calos with E>%d'%(E),100,0,20,deltaT)


        #calo hit position
        if event.isCaloHit:
            self.hist_svc.BookFillTH3(
                'caloHitPosition', 'first calo hit position', 100,0,200,100,0,200,100,0,200,
                event.caloHitX,event.caloHitY,event.caloHitZ)


        #wiggle plot for E>1.7 GeV
        self.fillWigglePlot(event)



    #bin size 149.2 ns
    #for 3k bins, maximum time reach to 3000*149.2 = 447.600 micro-senconds (about 7.5*muonLifeTime)
    #from tree->Draw, the time ranges approximatedly from 0 to 500e+3ns
    def fillWigglePlot(self,event):
        for time,E in event.time_energy_dict.items():
            #Time vs Energy
            ncalo = event.time_ncalo_dict[time]
            self.hist_svc.BookFillTH2('TE','Time vs Energy', 3000,0,447.6e3,1000,0,3000,time,E)
            self.hist_svc.BookFillTH1('T','N(t)',3000,0,447.6e3,time)
            self.hist_svc.BookFillTH1('T%d'%(ncalo),'N(t) : calo %d'%(ncalo),3000,0,447.6e3,time)

            energy_threshold = [1.4,1.5,1.6,1.7,1.8,1.9,2.0] # GeV
            for threshold in energy_threshold:
                if E>threshold:                    
                    self.hist_svc.BookFillTH1('Ethreshold%d_T'%(threshold*10),'N(t,E>%dGeV)'%(threshold),3000,0,447.6e3,time)
                    self.hist_svc.BookFillTH1('Ethreshold%d_T%d'%(threshold*10,ncalo),'N(t,E>%d GeV) : calorimeter %d'%(threshold,ncalo),3000,0,447.6e3,time)

    def fillHists(self,event):
        treeType = self.name_handler.GetSampleTag()
        eventType = 'omd' if event.posiTrackID == 4 else 'rmd'

        #set up event flags
        event_flags = {key:False for key in self.counter}
        event_flags['Total'] = True
        event_flags[treeType] = True                
        event_flags['%s_%s'%(treeType,eventType)] = True

        #fill cut flow
        for label,passed in event_flags.items():
            if passed:                
                self.hist_svc.BookFullCutHist('Cutflow','Cutflow',self.counter,label)

        #fill analysis hists
        self.fillAnalysisHists(event)
        self.name_handler.SetChannelTag(eventType)
        self.fillAnalysisHists(event)


def main():
    job = ReadVariables()
    job.WriteHistos()

if '__main__' == __name__:
    main()