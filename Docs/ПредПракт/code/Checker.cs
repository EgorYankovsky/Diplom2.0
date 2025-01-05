internal List<int> GetElem(double r, double z)
{
    int i = 0;
    for (; i < mesh2Dim.nodesR.Count - 1 && r >= 0.001; i++)
        if (mesh2Dim.nodesR[i] <= r && r <= mesh2Dim.nodesR[i + 1])
            break;
    int j = 0;
    for (; j < mesh2Dim.nodesZ.Count - 1; j++)
        if (mesh2Dim.nodesZ[j] <= z && z <= mesh2Dim.nodesZ[j + 1])
            break;
    return elemsArr[j * (mesh2Dim.nodesR.Count - 1) + i].Arr;
}