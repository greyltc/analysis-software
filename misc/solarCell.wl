(* ::Package:: *)

(* ::Input:: *)
(*SetOptions[$FrontEnd,"ClearEvaluationQueueOnKernelQuit"->False]*)
(*Quit[]*)


(* ::Input:: *)
(*(*Here's the charastic equation where y is current and x is voltage*)*)
(*chareqn=y==Iph-((x+y*Rs)/Rsh)-I0*(Exp[(x+y*Rs)/(n*Vth)]-1)*)


(* ::Input:: *)
(*(*What follows are the symbolic manipulations when we take x as the independant variable*)*)


(* ::Input:: *)
(*current = Solve[chareqn,y]*)


(* ::Input:: *)
(*Isc = Refine[y/.current,x==0]*)


(* ::Input:: *)
(*powerA=p==x*y/.current*)


(* ::Input:: *)
(*PPrimeA=D[x*y/.current,x]*)


(* ::Input:: *)
(*Vmpp = Solve[PPrimeA==0,x]*)


(* ::Input:: *)
(*(*What follows are the symbolic manipulations when we take y as the independant variable*)*)
(**)


(* ::Input:: *)
(*voltage = Solve[chareqn,x]*)


(* ::Input:: *)
(*Voc = Refine[x/.voltage,y==0]*)


(* ::Input:: *)
(*powerB=p==y*x/.voltage*)


(* ::Input:: *)
(*PPrimeB=D[y*x/.voltage,y]*)


(* ::Input:: *)
(*Impp = Solve[PPrimeB==0,y]*)



