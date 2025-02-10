**HEADING
ekate -> Abaqus interface
(c) 2016-2022, Hoang-Giang Bui, Ruhr-University Bochum
    2023-2024, Hoang-Giang Bui, Helmholtz-Zentrum Hereon
*if(strcmp(GenData(Export_Abaqus),"1")==0)
**Part, name=Ground
**NODE
*RealFormat "%10.10f"
*loop nodes
*nodesnum, *NodesCoord(1), *NodesCoord(2), *NodesCoord(3)
*end nodes
****
*loop elems
*tcl(WriteAbaqusInpElements *ElemsNum *ElemsType)
*end elems
****
*tcl(WriteAbaqusInpMaterials)
**End Part
*endif
