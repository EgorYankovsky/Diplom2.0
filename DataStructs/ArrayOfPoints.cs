namespace DataStructs;

public class ArrayOfPoints
{
    private List<Point> _list;

    public int Length => _list.Count;

    public Point? this[int i] => _list[i];

    public MyEnumerator GetEnumerator() => new(this);

    public class MyEnumerator
    {  
        int nIndex;  
        ArrayOfPoints collection;  
        public MyEnumerator(ArrayOfPoints coll)
        {  
            collection = coll;  
            nIndex = -1;  
        }  
    
        public bool MoveNext()
        {  
            nIndex++;  
            return nIndex < collection._list.Count;  
        }  
    
        public Point Current => collection._list[nIndex];
    }

    public void Append(Point p) => _list.Add(p);

    public ArrayOfPoints(string path)
    {
        _list = [];
        var data = File.ReadAllText(path).Split("\n");
        int pointsAmount = int.Parse(data[0]);
        for (int i = 0; i < pointsAmount; i++)
            _list.Add(new Point([.. data[i + 1].Split()]));
    }

    public ArrayOfPoints(int pointsAmount)
    {
        _list = new(pointsAmount);
    }
}