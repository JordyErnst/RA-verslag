# Verhoging van CD28 op T-helpercellen gemeten in Reumatoïde artritis patiënten. 
<p align="center">
  <img src="Assets/reuma-in-handen.jpg" alt="Reuma" width="2000" />
</p>

___

## Inleiding
Reumatoïde artritis (RA) is een auto-immuun ziekte. Het leidt tot een constante ontsteking van de gewrichtsvloeistof en afbraak van botten en kraakbeen, wat leidt tot de afbraak van de gewrichten [(Jahid et al., 2023)](Bronnen/(Jahid_et_al_2023).pdf).  De grootste oorzaak van reumatoïde artritis is genetisch, maar andere factoren zoals: geslacht, roken, BMI, alcoholconsumptie, dieet, en mondhygiëne verhogen ook de vatbaarheid voor RA [(Sparks, 2018)](Bronnen/(Sparks_2018).pdf).
Synoviumbiopten werden genomen van test en controle patienten, deze werden gesequencet en geannalyseerd binnen R met als doel het uitzoeken welke genen er hoger of lager tot expressie kwamen in personen met RA. Hiernaast werd er met behulp van een gene ontology analyses (GO-analyse) onderzocht welke meetabolsiche route het meest betrokken was bij RA.
In de onderstaande Flowschema is hoe dit experiment is uitgevoerd.


<p align="center">
  <img src="Assets/FlowSchema2.png" alt="Flowschema" width="700" />
</p>


## Methode
Tijdens dit onderzoek zijn 8 patiënten onderzocht. 4 ACPA (anti-ccp) negatief en 4 ACPA positief (RA). Synoviumbiopten van deze patiënten werden gesequenced en geanalyseerd met behulp van Rstudio. Binnen R werd er gebruik gemaakt van externe packages. Deze werden gedownload met `BiocManager (V 1.30.26)`.
`Rsubread (V2.20.0)` en `Rsamtools (V2.22.0)` werden gebruikt om een menselijk refentiegenoom `(Ensambl: GCA_000001405.29)` te indexeren. Op deze index werden de patienten sequenties gemapped. Hiermee werden reads van de patienten gekoppeld binnen aan locaties in het menselijkgenoom.
`Rsubread” (V2.20.0)` en `readr (V2.1.5)` werden gebruik om een count matrix te generenen. Een count matrix is een tabel wat aan toont hoeveel reads op elke gen zijn gemapt per patient.
Voor het laatste gedeelte werd er een(differentiële expressie analyse, KEGG-analyse, GO-analyse en volcano plot gegeneerd. Hiervoor waren de packages: `DESeq2 (V1.46.0)`, `KEGGREST (V1.46.0)`, `readr (V2.1.5)`, `dplyr (V1.1.4)`, `goseq (V1.58.0)`, `geneLenDataBase (V1.42.0)`, `org.Hs.eg.db (V3.20.0)`, `GO.db (V3.20.0)`, `ggplot2 (V3.5.2)`, `EnhancedVolcano (V1.24.0)` en `pathview (V1.46.0)` gebruikt. Aan de hand van deze analyses kon er onderscheid gemaakt worden tussen significante veranderingen binnen RA patiënten en kon er gemeten werden op welke metabole route RA het meest effect had.
Het volledige Rstudio code is [hier](Script/casus_Reuma_R_code.R) terug te vinden.


## Resultaten
Als eerste werd er een [volcano plot](Resultaten/volcano_plot.png) gegenereerd. Deze toont hoeveel genen niet significant (grijs), alleen Log2fc < 2 of >-2 hebben (groen) of ook een p-waarde van <0.05 hebben (rood). De GO-analyse genereerde een [top 10 grafiek](Resultaten/top_10_GO.png) van metabolische routes die het meest significant waren aangepast binnen RA. En als laatste werd de KEGG-analyse gebruikt om significante genen te visuliseren binnen de KEGG-pathway [“hsa04672”](Resultaten/hsa04672.png) (intenstinaal immuunnetwerk voor igA-productie). De pathway toont op gereguleerde genen aan (rood) en neer gereguleerde genen (groen). Hiermee is te zien dat er een duidelijk verschil is in genregulatie als gevolg van RA.

## Conclusie
Dit experiment werd uitgevoerd met als doel: uitzoeken welke genen er hoger of lager tot expressie komen in personen met reumatoïde artritis. In de KEGG-pathway is te zien dat de “T-cell receptor signaling pathway” opgereguleerd is. Specifiek het gen CD28+ op CD4-cellen (T-helpercellen), wat een grote rol speelt bij T-celactivatie. Uit eerdere onderzoek is gevonden dat een verhoogde klonale expansie van CD4+ CD28- T-cellen vaak wordt gevonden in patiënten met RA. [(Pawlik et al., 2003)](Bronnen/(Pawlik_et_al_2003).pdf)
Ook werd er een deelvraag gesteld: welke metabolische route het meest betrokken is bij patiënten met RA? Aan de hand van een GO-analyse is dit berekent. Binnen deze grafiek is een top 10 gegenereerd. Hierin is te zien dat het “immunoglobulinecomplex” het meest significant was aangepast binnen patiënten met RA.
