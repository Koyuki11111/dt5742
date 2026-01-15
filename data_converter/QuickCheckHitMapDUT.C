#include <filesystem>
#include <iostream>


using DAQ_KEY = pair<int,int>;
using STRIP_KEY = tuple<int,int,int>;

float amp_thr = 0.;

std::string to6digits(int n) {
  std::ostringstream oss;
  oss << std::setw(6) << std::setfill('0') << n;
  return oss.str();
}

bool CheckFileExist(string name){
  std::filesystem::path p = name;
  if (std::filesystem::exists(p)) return true;
  else                            return false; 
}

vector<float> GetLocalPos(int row, int col, bool ishorizon) {
  vector<float> pos(2);  
  if (ishorizon) {
    pos[0] = col+0.5;
    pos[1] = (row+1)*0.5;
  } else {
    pos[0] = (row+1)*0.5;
    pos[1] = col+0.5;    
  }

  return pos;
}

void QuickCheckHitMapDUT(int runnumber=1){

  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetStatBorderSize(0);
  gStyle->SetLegendBorderSize(0);

  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.1);
  
  string fname[2]; 
  TFile *input[2];
  TTree *ttree[2];

  for (int daq=0; daq<2; daq++) {
    fname[daq] = Form("/data/%s/daq0%d/output/root/waveforms_analyzed.root",to6digits(runnumber).c_str(),daq);
    if (!CheckFileExist(fname[daq])) {
      cout<<"[ERROR] The file "<<fname[daq]<<" doesn't exist."<<endl;
      cout<<"        Check run number or if you execute analyze_waveforms (stage-2) again!"<<endl;
      return ;
    } else {
      input[daq] = TFile::Open(fname[daq].c_str());
      ttree[daq] = (TTree*)input[daq]->Get("Analysis");
    }
  }
  
  if (ttree[0]->GetEntriesFast() != ttree[1]->GetEntriesFast()) {
    cout<<"[WORNING]  The number of events in daq00 and daq01 are not consistent."<<endl;
    cout<<"           Following analysis uses the one with fewer events."<<endl;
    cout<<"           daq00 is "<<ttree[0]->GetEntriesFast()<<" daq01 is "<<ttree[1]->GetEntriesFast()<<endl;
  } 
  
  int nEvt = std::min(ttree[0]->GetEntriesFast(),ttree[1]->GetEntriesFast());
  
  int __nChannels__=0;  
  std::vector<int>* __sensorID__=nullptr;
  std::vector<int>* __sensorCol__=nullptr;
  std::vector<int>* __sensorRow__=nullptr;
  std::vector<int>* __stripID__=nullptr;
  std::vector<bool>* __isHorizontal__=nullptr;
  std::vector<bool>* __hasSignal__=nullptr;
  std::vector<float>* __baseline__=nullptr;
  std::vector<float>* __rmsNoise__=nullptr;
  std::vector<float>* __noise1Point__=nullptr;
  std::vector<float>* __ampMinBefore__=nullptr;
  std::vector<float>* __ampMaxBefore__=nullptr;
  std::vector<float>* __ampMax__=nullptr;
  std::vector<float>* __charge__=nullptr;
  std::vector<float>* __signalOverNoise__=nullptr;
  std::vector<float>* __peakTime__=nullptr;
  std::vector<float>* __riseTime__=nullptr;
  std::vector<float>* __slewRate__=nullptr;
  std::vector<float>* __jitterRMS__=nullptr;
  float __jitterCFD__[2][16];
  
  int nChannels[2]={};  
  std::vector<int> sensorID[2];
  std::vector<int> sensorCol[2];
  std::vector<int> sensorRow[2];
  std::vector<int> stripID[2];
  std::vector<bool> isHorizontal[2];
  std::vector<bool> hasSignal[2];
  std::vector<float> baseline[2];
  std::vector<float> rmsNoise[2];
  std::vector<float> noise1Point[2];
  std::vector<float> ampMinBefore[2];
  std::vector<float> ampMaxBefore[2];
  std::vector<float> ampMax[2];
  std::vector<float> charge[2];
  std::vector<float> signalOverNoise[2];
  std::vector<float> peakTime[2];
  std::vector<float> riseTime[2];
  std::vector<float> slewRate[2];
  std::vector<float> jitterRMS[2];
  std::vector<float> jitterCFD[2];
  
  for (int daq=0; daq<2; daq++) {

    ttree[daq] -> SetBranchAddress("nChannels",&__nChannels__);
    ttree[daq] -> SetBranchAddress("sensorID",&__sensorID__);
    ttree[daq] -> SetBranchAddress("sensorCol",&__sensorCol__);
    ttree[daq] -> SetBranchAddress("sensorRow",&__sensorRow__);
    ttree[daq] -> SetBranchAddress("stripID",&__stripID__);
    ttree[daq] -> SetBranchAddress("isHorizontal",&__isHorizontal__);
    ttree[daq] -> SetBranchAddress("hasSignal",&__hasSignal__);
    ttree[daq] -> SetBranchAddress("baseline",&__baseline__);
    ttree[daq] -> SetBranchAddress("rmsNoise",&__rmsNoise__);
    ttree[daq] -> SetBranchAddress("noise1Point",&__noise1Point__);
    ttree[daq] -> SetBranchAddress("ampMinBefore",&__ampMinBefore__);
    ttree[daq] -> SetBranchAddress("ampMaxBefore",&__ampMaxBefore__);
    ttree[daq] -> SetBranchAddress("ampMax",&__ampMax__);
    ttree[daq] -> SetBranchAddress("charge",&__charge__);
    ttree[daq] -> SetBranchAddress("signalOverNoise",&__signalOverNoise__);
    ttree[daq] -> SetBranchAddress("peakTime",&__peakTime__);
    ttree[daq] -> SetBranchAddress("riseTime",&__riseTime__);
    ttree[daq] -> SetBranchAddress("slewRate",&__slewRate__);
    ttree[daq] -> SetBranchAddress("jitterRMS",&__jitterRMS__);
    for (int iCh=0; iCh<16; ++iCh) {
      string name= iCh<10 ? Form("0%d",iCh) : Form("%d",iCh);
      ttree[daq] -> SetBranchAddress(Form("ch%s_timeCFD_50pc",name.c_str()),&__jitterCFD__[daq][iCh]);
    }
  }

  
  TH2F * map_sensor[2];
  map_sensor[0] = new TH2F("map_sensor_vertical","",64,0,32,2,0,2);
  map_sensor[1] = new TH2F("map_sensor_horizontal","",2,0,2,64,0,32);

  TH2F* hist_amp_jitter = new TH2F("hist_amp_jitter","",3000,0,3000,500,0,0.1);
  hist_amp_jitter->SetTitle("; Amplitude [ADC]; Jitter [ns]");
  
  TH2F* hist_amp_risetime = new TH2F("hist_amp_risetime","",3000,0,3000,200,0,2);
  hist_amp_risetime->SetTitle("; Amplitude [ADC]; Rise-Time (10-90 pc) [ns]");

  TH2F* hist_peak_time = new TH2F("hist_peak_time","",3000,0,3000,200,0,2);
  hist_peak_time->SetTitle("; Amplitude [ADC]; Peak [ns]");

  
  TH2F* hist_hit_pos_sum = new TH2F("hist_hit_pos_sum","",400,-20,20,400,-20,20);
  TH2F* hist_hit_pos[4];
  for (auto idx=0; idx<4; idx++){
    hist_hit_pos[idx]= new TH2F(Form("hist_hit_pos_sensor%d",idx),"",400,-20,20,400,-20,20);
    hist_hit_pos[idx]->SetTitle(";Position [mm]; Position [mm]");
  }
  
  TH1F* hist_all_amp = new TH1F("hist_all_amp","",500,0,4000);
  hist_all_amp->SetTitle(";Amplitude [ADC];Entries");
  TH1F* hist_max_amp = new TH1F("hist_max_amp","",500,0,4000);
  hist_max_amp->SetTitle(";Amplitude [ADC];Entries");
  hist_max_amp->SetLineColor(kRed);
  
  TH1F* hist_diff_cfd_time[4][4];
  TH2F* hist_corr_xx[4][4];
  TH2F* hist_corr_yy[4][4];
  TH2F* hist_corr_xy[4][4];
  for (int sens1=0; sens1<4; sens1++) {
    for (int sens2=0; sens2<4; sens2++) {
      hist_diff_cfd_time[sens1][sens2] = new TH1F(Form("hist_diff_cfd_time_sensor%d_sensor%d",sens1,sens2),"",1000,0,100);
      hist_diff_cfd_time[sens1][sens2]->SetTitle(";#Delta CFD(50pc); Entries");

      hist_corr_xx[sens1][sens2] = new TH2F(Form("hist_corr_xx_sensor%d_sensor%d",sens1,sens2),"",400,-20,20,400,-20,20);
      hist_corr_xx[sens1][sens2]->SetTitle(Form(";x_%d [mm]; x_%d [mm]",sens1,sens2));      
      hist_corr_yy[sens1][sens2] = new TH2F(Form("hist_corr_yy_sensor%d_sensor%d",sens1,sens2),"",400,-20,20,400,-20,20);
      hist_corr_yy[sens1][sens2]->SetTitle(Form(";y_%d [mm]; y_%d [mm]",sens1,sens2));      
      hist_corr_xy[sens1][sens2] = new TH2F(Form("hist_corr_xy_sensor%d_sensor%d",sens1,sens2),"",400,-20,20,400,-20,20);
      hist_corr_xy[sens1][sens2]->SetTitle(Form(";x_%d [mm]; y_%d [mm]",sens1,sens2));      
    }
  }


  
  
  int numV=0;
  int numH=0;
  bool isChecked=false;
  
  for (int iEvt=0; iEvt<nEvt; ++iEvt) {
    
    set<int> unique_sensorIDs;
    
    for (int daq=0; daq<2; daq++) {
      jitterCFD[daq].clear();
      ttree[daq] -> GetEntry(iEvt);      
      nChannels[daq] = __nChannels__;
      sensorID[daq] = *__sensorID__;
      sensorCol[daq] = *__sensorCol__;
      sensorRow[daq] = *__sensorRow__;
      isHorizontal[daq] = *__isHorizontal__;
      baseline[daq] = *__baseline__;
      rmsNoise[daq] = *__rmsNoise__;
      ampMax[daq] = *__ampMax__;
      signalOverNoise[daq] = *__signalOverNoise__;
      peakTime[daq] = *__peakTime__;
      riseTime[daq] = *__riseTime__;
      slewRate[daq] = *__slewRate__;
      jitterRMS[daq] = *__jitterRMS__;

      for (int iCh=0; iCh<16; ++iCh) {
	jitterCFD[daq].push_back(__jitterCFD__[daq][iCh]);
      }

      for (auto ids : sensorID[daq]) {
	unique_sensorIDs.insert(ids);
      }
    }
    
    map<DAQ_KEY,STRIP_KEY> map_daq2strip;
    map<STRIP_KEY,DAQ_KEY> map_strip2daq;
    
    map<int, bool> map_ishorizontal_sensor;    
    map<int,map<pair<int,int>,float>> hitMap;
    map<int,DAQ_KEY> high_amp_daqid;
    map<int, float> high_amp;
    high_amp[0] = 0.0;
    high_amp[1] = 0.0;
    high_amp[2] = 0.0;
    high_amp[3] = 0.0;
    
    for (int daq=0; daq<2; daq++) {
      for (int iCh=0; iCh<nChannels[daq]; ++iCh) {
	int sens_id = sensorID[daq].at(iCh);	
	int row = sensorRow[daq].at(iCh);
	int col = sensorCol[daq].at(iCh);
	float amp = ampMax[daq].at(iCh);	
	map_ishorizontal_sensor[sens_id] = isHorizontal[daq].at(iCh);
	hitMap[sens_id][{row,col}]       = amp;	
	map_daq2strip[{daq,iCh}]         = {sens_id,row,col};
	map_strip2daq[{sens_id,row,col}] = {daq,iCh};
	
	if (amp>high_amp[sens_id]) {
	  high_amp[sens_id] = amp;	  
	  DAQ_KEY daq_key = map_strip2daq[{sens_id,row,col}];
	  high_amp_daqid[sens_id] = daq_key;
	}

      }
    }

    map<int, pair<float,float>> map_sensor_hitpos;

    float event_hitpos[2]={};
    int nH=0;
    int nV=0;

    float max_cfd_time[4]={}; 
    float meas_x[4]={};
    float meas_y[4]={};
    
    for (auto sens_id : unique_sensorIDs) {
      int max_daq=high_amp_daqid[sens_id].first;
      int max_ch =high_amp_daqid[sens_id].second;
      float max_amp = high_amp[sens_id];

      max_cfd_time[sens_id] = jitterCFD[max_daq].at(max_ch); 
      
      hist_amp_jitter->Fill(max_amp,jitterRMS[max_daq].at(max_ch));  
      hist_amp_risetime->Fill(max_amp,riseTime[max_daq].at(max_ch));  
      hist_max_amp->Fill(max_amp);
      hist_peak_time->Fill(max_amp,peakTime[max_daq].at(max_ch));
      
      bool isH = isHorizontal[max_daq].at(max_ch);
      
      const auto map_amp = hitMap[sens_id];
      
      float mean_row = 0;
      float mean_col = sensorCol[max_daq].at(max_ch) == 0 ? 5 : 15;
      float sum_amp  = 0;
      
      for (const auto &[rowcol, amp] : map_amp) {
	int row = rowcol.first;
	int col = rowcol.second;
	
	if (col == sensorCol[max_daq].at(max_ch)) {
	  vector<float> pos = GetLocalPos(row,col,isHorizontal[max_daq].at(max_ch));
	  hist_all_amp->Fill(amp);
	  
	  if (amp<amp_thr) continue;
	  if (isH){	    
	    mean_row += pos[1]*amp;
	    sum_amp  += amp;
	  } else {
	    mean_row += pos[0]*amp;
	    sum_amp  += amp;	    
	  }	  
	}
      }
      
      mean_row /= sum_amp;
      
      if (isH) {	
	map_sensor_hitpos[sens_id]={mean_col,mean_row};	
	hist_hit_pos[sens_id]->Fill(mean_col-10,mean_row-16);
	meas_x[sens_id] = mean_col-10;
	meas_y[sens_id] = mean_row-16;
	event_hitpos[1] += mean_row;
	nH++;
      } else {
	map_sensor_hitpos[sens_id]={mean_row,mean_col};
	hist_hit_pos[sens_id]->Fill(mean_row-16,mean_col-10);
	meas_x[sens_id] = mean_row-16;
	meas_y[sens_id] = mean_col-10;
	event_hitpos[0] += mean_row;
	nV++;
      }            
    }

    for (auto sens_id1 : unique_sensorIDs){
      for (auto sens_id2 : unique_sensorIDs){
	hist_diff_cfd_time[sens_id1][sens_id2]->Fill(max_cfd_time[sens_id1]-max_cfd_time[sens_id2]);
      }
    }

    for (auto sens_id1 : unique_sensorIDs){
      for (auto sens_id2 : unique_sensorIDs){
	hist_corr_xx[sens_id1][sens_id2]->Fill(meas_x[sens_id1],meas_x[sens_id2]);	
	hist_corr_yy[sens_id1][sens_id2]->Fill(meas_y[sens_id1],meas_y[sens_id2]);	
	hist_corr_xy[sens_id1][sens_id2]->Fill(meas_x[sens_id1],meas_y[sens_id2]);
      }
    }
    
    if (!isChecked) {
      numV = nV;
      numH = nH;
      isChecked=true;
    }
    
    if (nV>0) {
      event_hitpos[0]/=nV;
      event_hitpos[0] -= 16;
    }
    else {
      event_hitpos[0] = 5;
      event_hitpos[0] -= 10;
    }    
    if (nH>0) {
      event_hitpos[1]/=nH;
      event_hitpos[1] -= 16;
    } else {
      event_hitpos[1] = 5;
      event_hitpos[1] -= 10;
    }
    
    hist_hit_pos_sum->Fill(event_hitpos[0],event_hitpos[1]);
  }
  
  TLine* lineH = new TLine(0,-20,0,20);
  lineH->SetLineWidth(2);
  lineH->SetLineStyle(2);
  lineH->SetLineColor(kBlack);

  TLine* lineV = new TLine(-20,0,20,0);
  lineV->SetLineWidth(2);
  lineV->SetLineStyle(2);
  lineV->SetLineColor(kBlack);

  TBox *boxH = new TBox(-10, -16, 10, 16);
  boxH->SetFillStyle(0);
  boxH->SetLineWidth(2);
  boxH->SetLineColor(kBlack);
  
  TBox *boxV = new TBox(-16, -10, 16, 10);
  boxV->SetFillStyle(0);
  boxV->SetLineWidth(2);
  boxV->SetLineColor(kBlack);

  TLine* lineHitMeanV = new TLine(hist_hit_pos_sum->GetMean(1),-20,hist_hit_pos_sum->GetMean(1),20);
  lineHitMeanV->SetLineWidth(2);
  lineHitMeanV->SetLineStyle(1);
  lineHitMeanV->SetLineColor(kRed);

  TLine* lineHitMeanH = new TLine(-20,hist_hit_pos_sum->GetMean(2),20,hist_hit_pos_sum->GetMean(2));
  lineHitMeanH->SetLineWidth(2);
  lineHitMeanH->SetLineStyle(1);
  lineHitMeanH->SetLineColor(kRed);

  TText *txtV = new TText(-18.53812,17.35126, Form("Mean-X = %.4f [mm]",hist_hit_pos_sum->GetMean(1)));
  txtV->SetTextSize(0.04);
  txtV->SetTextColor(kRed);

  TText *txtH = new TText(1.801594,-1.902521, Form("Mean-Y = %.4f [mm]",hist_hit_pos_sum->GetMean(2)));
  txtH->SetTextSize(0.04);
  txtH->SetTextColor(kRed);

  if (numH==0) {
    hist_hit_pos[0]->SetTitle("Weighted Mean Position (Only Vertical Sensor)");
  } else if (numV==0) {  
    hist_hit_pos[0]->SetTitle("Weighted Mean Position (Only Hrizontal Sensor)");
  } else {
    hist_hit_pos[0]->SetTitle(Form("Weighted Mean Position (%d x Hrizontal and %d x Vertical Sensor)",numH,numV));
  }
  
  TCanvas* c1 = new TCanvas("c1","",900,900);
  c1->cd(1);
  hist_hit_pos[0]->Draw("col");
  hist_hit_pos[1]->Draw("col same");
  hist_hit_pos_sum->Draw("col same");
  if (numH>0) boxH->Draw("same");
  if (numV) boxV->Draw("same");
  lineH->Draw("same");
  lineV->Draw("same");
  lineHitMeanH->Draw("same");
  lineHitMeanV->Draw("same");
  txtV->Draw("same");
  txtH->Draw("same");

  TCanvas* c2 = new TCanvas("c2","",900,900);
  c2->Divide(2,2);
  c2->cd(1);
  hist_amp_jitter->Draw("col");
  c2->cd(2);
  hist_amp_risetime->Draw("colz");
  c2->cd(3);
  gPad->SetLogy(1);
  hist_all_amp->Draw("");
  hist_max_amp->Draw("same");
  c2->cd(4);
  hist_diff_cfd_time[0][1]->Draw("");

  TCanvas* c3 = new TCanvas("c3","",1200,300);
  c3->Divide(4,1);
  c3->cd(1);
  hist_corr_xx[0][1]->Draw("");
  c3->cd(2);
  hist_corr_yy[0][1]->Draw("");
  c3->cd(3);
  hist_corr_xy[0][1]->Draw("");
  c3->cd(4);
  hist_corr_xy[1][0]->Draw("");

  
  cout<<"The number of Horizontal Sensor = "<<numH<<endl;
  cout<<"The number of Vertical Sensor = "<<numV<<endl;
}
