run_5:
Germinal VHH design against GPC3 hotspots:
Used structure correct_uniprot_gpc3_59-477_hotspot2.pdb
The goal is to design binders to the outer surface of GPC3 rather than the internal surface that run_1 has produced to.
Hotspots: 301,467,471,393,400, Renumbered: A22,A73,A77,A42,A49
Parameters changed for VHH and ScFv runs -- mostly from IL13 and PDL1 runs.
Four independent runs launched, one on each GPU:{0,3}.
No successful designs. But several reach the AbMPNN stage. See: ./results/trajectories/ .
However the structures in the trajectories directory are not great (weird binding modes for an antibody).
Also Germinal is not respecting the target template. Thus abandon.
