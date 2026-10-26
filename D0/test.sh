root -l -b -q -e '
TFile f("/Users/ml2024/SilasThesis/Masters-Thesis/D0/data/002/tof-qa/AnalysisResults.root");
if(f.IsZombie())gSystem->Exit(1);
gStyle->SetOptFit(111);
TCanvas canvas("tofcheck","TOF response",1200,1200);
canvas.Divide(2,4);
double edges[]={0.5,1.0,2.0,3.0,5.0};
int column=0;
for(auto species:{"Pi","Ka"}){
  ++column;
  auto h2=f.Get<TH2>(Form("pid-tof-qa-mc/nsigmaMCgoodmatch/%s/%s",species,species));
  if(!h2){std::cout<<species<<": HISTOGRAM_MISSING\n";continue;}
  for(int i=0;i<4;++i){
    canvas.cd(2*i+column);
    int first=h2->GetXaxis()->FindFixBin(edges[i]+1e-6);
    int last=h2->GetXaxis()->FindFixBin(edges[i+1]-1e-6);
    auto h=h2->ProjectionY(Form("core_%s_%d",species,i),first,last);
    h->SetTitle(Form("%s: %.1f < pT < %.1f GeV/c",species,edges[i],edges[i+1]));
    h->GetXaxis()->SetRangeUser(-5,5);
    h->Draw();
    double core=h->Integral(h->FindBin(-2),h->FindBin(2));
    if(core>=100){
      int status=h->Fit("gaus","Q","",-2,2);
      if(status)std::cout<<species<<" slice "<<i<<": FIT_FAILED\n";
    }else{
      std::cout<<species<<" slice "<<i<<": INSUFFICIENT_CORE_STATISTICS\n";
    }
  }
}
canvas.SaveAs("/Users/ml2024/SilasThesis/Masters-Thesis/D0/data/002/tof-qa/TOF-core-check.pdf");
'