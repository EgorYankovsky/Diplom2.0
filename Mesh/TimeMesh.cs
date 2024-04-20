namespace Grid;

public class TimeMesh(double[] timeLayers)
{
    private double[] _timeLayers = timeLayers;

    public double this[int i] => _timeLayers[i]; 

    public int Count => _timeLayers.Length;
}