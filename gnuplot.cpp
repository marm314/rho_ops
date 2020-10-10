#include"gnuplot.h"

GNUPLOT::GNUPLOT(){path="";}
GNUPLOT::~GNUPLOT()
{}
void GNUPLOT::plot3D(string data,string gnuplot_commands,string operation,string title_op,string *extracommands,
int extra_commands,char var1,char var2,double var1_r[2],double var2_r[2])
{
 int i;
 quotation="\"";
 ofstream commands;
 //PM3D
 commands.open((gnuplot_commands+".plot3").c_str());
 commands<<"set term postscrip eps color enhanced "+quotation+"Helvetica"+quotation+" 13"<<endl;
 commands<<"set output "+quotation+operation+".eps"+quotation<<endl;
 commands<<"set title "+quotation+title_op+" vs ("<<var1<<","<<var2<<",0)"+quotation+" font "+quotation+"Helvetica,18"+quotation<<endl;
 commands<<"set xlabel "+quotation<<var1<<quotation+" font "+quotation+"Helvetica,18"+quotation<<endl;
 commands<<"set ylabel "+quotation<<var2<<quotation+" font "+quotation+"Helvetica,18"+quotation<<endl;
 commands<<"set zlabel "+quotation+title_op+quotation+" font "+quotation+"Helvetica,18"+quotation+"offset 2,11"<<endl;
 commands<<"set xr ["<<var1_r[0]<<":"<<var1_r[1]<<"]"<<endl;
 commands<<"set yr ["<<var2_r[0]<<":"<<var2_r[1]<<"]"<<endl;
 commands<<"set format z "+quotation+"%.2e"+quotation<<endl;
 commands<<"set grid"<<endl;
 commands<<"set contour base"<<endl;
 commands<<"set ticslevel 0"<<endl;
 commands<<"set key spacing 1"<<endl;
 commands<<"set isosamples 40"<<endl;
 commands<<"set cntrparam levels 10"<<endl;
 commands<<"unset clabel"<<endl;
 commands<<"set style line 1 lt rgb "+quotation+"black"+quotation+"lw 1.5 pt 4"<<endl;
 if(extra_commands!=0)
 {
  for(i=0;i<extra_commands;i++)
  {
   commands<<extracommands[i]<<endl;
  }
 }
 commands<<"splot "+quotation+data+quotation+" u 1:2:3 w pm3d ls 1 title "+quotation+" "+quotation<<endl;
 commands.close();
 if(path!="")
 {
  system((path+" <"+gnuplot_commands+".plot3").c_str());
 }
 else
 {
  system(("gnuplot <"+gnuplot_commands+".plot3").c_str());
 }
 //PM3D MAP
 commands.open((gnuplot_commands+"_map.plot3").c_str());
 commands<<"set term postscrip eps color enhanced "+quotation+"Helvetica"+quotation+" 13"<<endl;
 commands<<"set output "+quotation+operation+"_map.eps"+quotation<<endl;
 commands<<"set title "+quotation+title_op+" vs ("<<var1<<","<<var2<<",0)"+quotation+" font "+quotation+"Helvetica,20"+quotation<<endl;
 commands<<"set xlabel "+quotation<<var1<<quotation+" font "+quotation+"Helvetica,20"+quotation<<endl;
 commands<<"set label 1 "+quotation<<var2<<quotation+" font "+quotation+"Helvetica,20"+quotation+" at graph -0.1, graph 0.5"<<endl;
 commands<<"set zlabel "+quotation+title_op+quotation+" font "+quotation+"Helvetica,20"+quotation<<endl;
 commands<<"set xr ["<<var1_r[0]<<":"<<var1_r[1]<<"]"<<endl;
 commands<<"set yr ["<<var2_r[0]<<":"<<var2_r[1]<<"]"<<endl;
 commands<<"set grid"<<endl;
 commands<<"set contour base"<<endl;
 commands<<"set pm3d map"<<endl;
 commands<<"set key spacing 1"<<endl;
 commands<<"set palette rgbformulae 33,13,10"<<endl;
 commands<<"set cntrparam levels 10"<<endl;
 commands<<"unset clabel"<<endl;
 commands<<"set style line 1 lt rgb "+quotation+"black"+quotation+"lw 1.5 pt 4"<<endl;
 if(extra_commands!=0)
 {
  for(i=0;i<extra_commands;i++)
  {
   commands<<extracommands[i]<<endl;
  }
 }
 commands<<"splot "+quotation+data+quotation+" u 1:2:3 w pm3d ls 1 title "+quotation+" "+quotation<<endl;
 commands.close();
 if(path!="")
 {
  system((path+" <"+gnuplot_commands+"_map.plot3").c_str());
 }
 else
 {
  system(("gnuplot <"+gnuplot_commands+"_map.plot3").c_str());
 }
}

void GNUPLOT::plot2D(string data,string gnuplot_commands,string operation,string title_op,string *extracommands,
int extra_commands,char var1,double var1_r[2])
{
 int i;
 quotation="\"";
 ofstream commands;
 //2D Plot
 commands.open((gnuplot_commands+".plot2").c_str());
 commands<<"set term postscrip eps color enhanced "+quotation+"Helvetica"+quotation+" 13"<<endl;
 commands<<"set output "+quotation+operation+".eps"+quotation<<endl;
 commands<<"set title "+quotation+title_op+" vs "<<var1<<" "<<quotation+" font "+quotation+"Helvetica,18"+quotation<<endl;
 commands<<"set xlabel "+quotation<<var1<<quotation+" font "+quotation+"Helvetica,18"+quotation<<endl;
 commands<<"set label 1 "+quotation+title_op+quotation+" font "+quotation+"Helvetica,18"+quotation+" at graph -0.05, graph 1.06"<<endl;
 commands<<"set xr ["<<var1_r[0]<<":"<<var1_r[1]<<"]"<<endl;
 commands<<"set grid"<<endl;
 commands<<"set format y "+quotation+"%.2e"+quotation<<endl;
 if(extra_commands!=0)
 {
  for(i=0;i<extra_commands;i++)
  {
   commands<<extracommands[i]<<endl;
  }
 }
 commands<<"plot "+quotation+data+quotation+" u 1:2 with linespoints title "+quotation+title_op+quotation<<endl;
 commands.close();
 if(path!="")
 {
  system((path+" <"+gnuplot_commands+".plot2").c_str());
 }
 else
 {
  system(("gnuplot <"+gnuplot_commands+".plot2").c_str());
 }
 commands.close();
}
