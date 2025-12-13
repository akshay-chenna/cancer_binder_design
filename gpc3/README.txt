1. structures: starting structures for binder design

2. run_1: bindcraft minibinder design

3. run_2:
Minibinder design runs with hotspots and trimmed structure.
/paperspace/Desktop/cancer/gpc3/run_2
Used structure correct_uniprot_gpc3_59-477_hotspot1.pdb
The goal is to design binders to the outer surface of GPC3 rather than the internal surface that run_1 has produced to.
Hotspots: 117, 125, 214, 311, 313, 314, 324, 333
Hotspots not respected. Killed.

4. run_3:
Minibinder design runs with hotspots and trimmed structure.
/paperspace/Desktop/cancer/gpc3/run_3
Used structure correct_uniprot_gpc3_59-477_hotspot2.pdb
The goal is to design binders to the outer surface of GPC3 rather than the internal surface that run_1 has produced to.
Hotspots: 301,467,471,393,400
Hotspots not respected. Killed.

5. run_4:
Germinal VHH design against GPC3 hotspots:
Used structure correct_uniprot_gpc3_59-477_hotspot1.pdb
The goal is to design binders to the outer surface of GPC3 rather than the internal surface that run_1 has produced to.
Hotspots: 117, 125, 214, 311, 313, 314, 324, 333
No progression to AbMPNN stage. Killed.

6. run_5:
Germinal VHH design against GPC3 hotspots:
Used structure correct_uniprot_gpc3_59-477_hotspot2.pdb
The goal is to design binders to the outer surface of GPC3 rather than the internal surface that run_1 has produced to.
Hotspots: 301,467,471,393,400
Hotspots not respected. Killed.

7. run_6:
Germinal ScFv design against GPC3 hotspots:
Used structure correct_uniprot_gpc3_59-477_hotspot2.pdb
The goal is to design binders to the outer surface of GPC3 rather than the internal surface that run_1 has produced to.
Hotspots: 301,467,471,393,400
Insufficient GPU memory.

8. run_7:
mBER VHH design against GPC3 hotspots:
Used structure correct_uniprot_gpc3_59-477_hotspot2.pdb
The goal is to design binders to the outer surface of GPC3 rather than the internal surface that run_1 has produced to.
Hotspots: 301,467,471,393,400

9. run_8:
BoltzGen minibinder design against GPC3 substructure 1 with hotspots. Binders to the outer surface of GPC3 rather than the interal surface.
Structure: Full structure input but visible: correct_uniprot_gpc3_59-477_renum.pdb (renumbered from 1) 
Hotspots: 67,71,147,153,156,275 (renumbered). Hotspots from SAP score. See 7zaw_sap_annotated_h1.pse
Designs with and without cysteins

10. run_9:
BoltzGen VHH design against GPC3 substructure 1 with hotspots. Binders to the outer surface of GPC3 rather than the interal surface.
Structure: Full structure input but visible: correct_uniprot_gpc3_59-477_renum.pdb (renumbered from 1) 
Hotspots: 67,71,147,153,156,275 (renumbered). Hotspots from SAP score. See 7zaw_sap_annotated_h1.pse

10. run_10:
BoltzGen minibinder design against GPC3 substructure 2 with hotspots. Binders to the outer surface of GPC3 rather than the interal surface.
Structure: Full structure input but visible: correct_uniprot_gpc3_59-477_renum.pdb (renumbered from 1)
Hotspots: 232,242,243,404,405 (renumbered). Hotspots from SAP score. See 7zaw_sap_annotated_h2.pse
Designs with and without cysteins.

11. run_11:
BoltzGen VHH design against GPC3 substructure 2 with hotspots. Binders to the outer surface of GPC3 rather than the interal surface.
Structure: Full structure input but visible: correct_uniprot_gpc3_59-477_renum.pdb (renumbered from 1)
Hotspots: 232,242,243,404,405 (renumbered). Hotspots from SAP score. See 7zaw_sap_annotated_h2.pse

12. run_12:
Modal replica of run_8. No restrictiions on cysteins. 20000 designs.
BoltzGen minibinder design against GPC3 substructure 1 with hotspots. Binders to the outer surface of GPC3 rather than the interal surface.
Structure: Full structure input but visible: correct_uniprot_gpc3_59-477_renum.pdb (renumbered from 1) 
Hotspots: 67,71,147,153,156,275 (renumbered). Hotspots from SAP score. See 7zaw_sap_annotated_h1.pse

13. run_13:
Modal replica of run_9. 20000 designs.
BoltzGen VHH design against GPC3 substructure 1 with hotspots. Binders to the outer surface of GPC3 rather than the interal surface.
Structure: Full structure input but visible: correct_uniprot_gpc3_59-477_renum.pdb (renumbered from 1) 
Hotspots: 67,71,147,153,156,275 (renumbered). Hotspots from SAP score. See 7zaw_sap_annotated_h1.pse

14. run_14:
Modal replica of run_10. No restriction on cysteins. 7500 designs.
BoltzGen minibinder design against GPC3 substructure 2 with hotspots. Binders to the outer surface of GPC3 rather than the interal surface.
Structure: Full structure input but visible: correct_uniprot_gpc3_59-477_renum.pdb (renumbered from 1)
Hotspots: 232,242,243,404,405 (renumbered). Hotspots from SAP score. See 7zaw_sap_annotated_h2.pse
Designs with and without cysteins.

15. run_15: # Not ordered
Modal replica of run_11: 7500 designs. # Not ordered
BoltzGen VHH design against GPC3 substructure 2 with hotspots. Binders to the outer surface of GPC3 rather than the interal surface.
Structure: Full structure input but visible: correct_uniprot_gpc3_59-477_renum.pdb (renumbered from 1)
Hotspots: 232,242,243,404,405 (renumbered). Hotspots from SAP score. See 7zaw_sap_annotated_h2.pse

16: run_16:
Protein Hunter minibinder run without cif/pdb templates.
Hotspot type 1: Same as in run_12 (from SAP score). 67,71,147,153,156,275 (renumbered)
PDB to use: correct_uniprot_gpc3_59-477_renum.pdb (renumbered)
Does not use templates.

17. run_17:
Protein Hunter minibinder run without cif/pdb templates.
Hotspot type 2: Same as in run_14 (from SAP score). 232,242,243,404,405 (renumbered)
PDB to use: correct_uniprot_gpc3_59-477_renum.pdb (renumbered)

18. run_18:
BEST IS TO IGNORE THIS RUN. PREFER RUN_21 INSTEAD.
RFD3 minibinder using correct_uniprot_gpc3_59-477_renum.pdb (renumbered from 1).
Hotspot type 1: Same as in run_12 (from SAP score). 67,71,147,153,156,275 (renumbered)
5000 designs

19. run_19:
BEST IS TO IGNORE THIS RUN. PREFER RUN_22 INSTEAD.
RFD3 minibinder using correct_uniprot_gpc3_59-477_renum.pdb (renumbered from 1).
Hotspot type 2: Same as in run_14 (from SAP score). 232,242,243,404,405 (renumbered)
2000 designs

20. run_20:
RFAntibody using YP7 scaffold.

21. run_21:
RFD3 minibinder using correct_uniprot_gpc3_59-477_renum.pdb (renumbered from 1).
Hotspot type 1: Same as in run_12 (from SAP score). 67,71,147,153,156,275 (renumbered)
10000 designs. MPNN temperature from RFD1 SI. SolubleMPNN weights. n_mpnn=6, larger minibinders=80-200

22. run_22:
RFD3 minibinder using correct_uniprot_gpc3_59-477_renum.pdb (renumbered from 1).
Hotspot type 2: Same as in run_14 (from SAP score). 232,242,243,404,405 (renumbered)
10000 designs. MPNN temperature from RFD1 SI. SolubleMPNN weights. n_mpnn=6, larger minibinders=80-200
