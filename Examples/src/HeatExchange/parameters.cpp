#include <iostream>
#include "parameters.hpp"
std::ostream & operator << (std::ostream & out,const parameters & p)
{
  out<<"PARAMETER VALUES:"<<"\n";
  out<<"L= "<<p.L<<"\n";
  out<<"a1= "<<p.a1<<"\n";
  out<<"a2= "<<p.a2<<"\n";
  out<<"To= "<<p.To<<"\n";
  out<<"Te= "<<p.Te<<"\n";
  out<<"k= "<<p.k<<"\n";
  out<<"hc= "<<p.hc<<"\n";
  out<<"M= "<<p.M<<"\n";
  out<<"nome del file dei risultati: "<<p.newname<<"\n";
  out<<"formato dell'output: "<<p.output<<"\n\n";
  return out;
}
