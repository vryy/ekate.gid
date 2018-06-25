**HEADING
ekate -> Abaqus interface
(c) 2016, 2017, 2018 Hoang-Giang Bui, Ruhr-University Bochum
**NODE
*RealFormat "%10.5f"
*loop nodes
*nodesnum, *NodesCoord(1), *NodesCoord(2), *NodesCoord(3)
*end nodes
****
*loop elems
*tcl(WriteAbaqusInpElements *ElemsNum *ElemsType)
*end elems
****
*tcl(WriteAbaqusInpMaterials)

