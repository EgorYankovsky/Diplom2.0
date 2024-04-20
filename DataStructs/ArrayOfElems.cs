namespace DataStructs;

public class ArrayOfElems
{
    public List<List<int>> Arr;
    
    public List<double> mui;

    public List<double> sigmai;

    public int Length { get => Arr.Count; set {} }

    public MyEnumerator GetEnumerator() => new(this);

    public List<int> this[int i] => Arr[i];

    public class MyEnumerator
    {  
        int nIndex;  
        ArrayOfElems collection;  
        public MyEnumerator(ArrayOfElems coll)
        {  
            collection = coll;  
            nIndex = -1;  
        }  
    
        public bool MoveNext()
        {  
            nIndex++;  
            return nIndex < collection.Arr.Count;  
        }  
    
        public List<int> Current => collection.Arr[nIndex];
    }

    public ArrayOfElems(string path)
    {
        Arr = [];
        mui = [];
        sigmai = [];

        var data = File.ReadAllText(path).Split("\n");
        int length = int.Parse(data[0]);
        for (int i = 0; i < length; i++)
        {
            var info = data[1 + i].Split(" ");
            Arr.Add([int.Parse(info[0]), int.Parse(info[1]), int.Parse(info[2]), int.Parse(info[3])]);
            mui.Add(double.Parse(info[4]));
            sigmai.Add(double.Parse(info[5]));
        }
    }


    public void Add(List<int> elem)
    {
        Arr.Add(elem);
    }

    public ArrayOfElems(int elemsAmount)
    {
        Arr = new(elemsAmount);
        mui = new();
        sigmai = new();
    }
}