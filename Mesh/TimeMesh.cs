using System.Collections;

namespace Grid;

public class TimeMesh(double[] timeLayers) : IEnumerable
{
    private double[] _timeLayers = timeLayers;

    public double this[int i] => _timeLayers[i]; 

    public int Count => _timeLayers.Length;

    public IEnumerator GetEnumerator() => _timeLayers.GetEnumerator();
}