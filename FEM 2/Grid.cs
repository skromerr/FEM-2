using System.Collections;
using System.Collections.Specialized;
using System.Drawing;

namespace UMFCourseProject;

public struct Material
{
   public double Mu { get; init; }
   public double J { get; init; }
   public string Name { get; init; }

   public Material(double mu, double j, string name)
   {
      Mu = mu;
      J = j;
      Name = name;
   }
}

public class Grid
{
    public List<Point2D> Nodes { get; set; }
    public List<int> FirstConditionNodes { get; init; }
    public int[][] Elements { get; init; }
    public Material[] Materials { get => _materials; }

    private Material[] _materials;

    private double _maxX, _maxY, _minX;

    public Grid(string spaceGridPath, string materialsPath, bool isQuadraticBasis = false)
    {
        string coordsPath = spaceGridPath + "coords.txt";
        string elemsPath = spaceGridPath + "elemsTr.txt";
        string firstCondPath = spaceGridPath + "firstConds.txt";
        string tokuPath = materialsPath + "toku";
        string muPath = materialsPath + "mu";

        setMaterials(muPath, tokuPath);
        string[] data;

        using (StreamReader sr = new(coordsPath))
        {


            data = sr.ReadLine()!.Split().ToArray();
            var nodesCount = Convert.ToInt32(data[0]);
            Nodes = new List<Point2D>();
            for (int i = 0; i < nodesCount; i++)
            {
                data = sr.ReadLine()!.Split(" ").ToArray();
                Nodes.Add(new Point2D(Convert.ToDouble(data[0]), Convert.ToDouble(data[1])));
            }

        }

        using (StreamReader sr = new(elemsPath))
        {

            data = sr.ReadLine()!.Split().ToArray();

            if (!isQuadraticBasis)
                Elements = new int[Convert.ToInt32(data[0])].Select(_ => new int[4]).ToArray(); //линейный
            else 
                Elements = new int[Convert.ToInt32(data[0])].Select(_ => new int[7]).ToArray(); // квадратичный
            
            for (int i = 0; i < Elements.Length; i++)
            {
                data = sr.ReadLine()!.Split(" ").ToArray();
                Elements[i][0] = Convert.ToInt32(data[0]) - 1;
                Elements[i][1] = Convert.ToInt32(data[1]) - 1;
                Elements[i][2] = Convert.ToInt32(data[2]) - 1;
                Elements[i][^1] = Convert.ToInt32(data[3]) - 1;
            }

        }

        using (StreamReader sr = new(firstCondPath))
        {

            data = sr.ReadToEnd()!.Split("\n").Where(str => str != "").ToArray();
            FirstConditionNodes = new List<int>();
            int nodeNum = Convert.ToInt32(data[1]) - 1;
            _maxX = Nodes[nodeNum].X;
            _maxY = Nodes[nodeNum].Y;
            _minX = Nodes[nodeNum].X;
            FirstConditionNodes.Add(nodeNum);


            for (int i = 2; i < data.Length; i++)
            {
                nodeNum = Convert.ToInt32(data[i]) - 1;
                FirstConditionNodes.Add(nodeNum);

                if (_maxX < Nodes[nodeNum].X) _maxX = Nodes[nodeNum].X;
                if (_maxY < Nodes[nodeNum].Y) _maxY = Nodes[nodeNum].Y;
                if (_minX > Nodes[nodeNum].X) _minX = Nodes[nodeNum].X;
            }

            

            if (FirstConditionNodes.Count != Convert.ToInt32(data[0])) { throw new NotImplementedException(); }

            
        }

        if (isQuadraticBasis)
        {
            //int count1 = FirstConditionNodes.Count;
            //Console.WriteLine($"Количество узлов под 1 краевым до добавления узлов от квадратичного базиса - {FirstConditionNodes.Count}");

            NumberNodes(); // для квадратичного базиса

            //Console.WriteLine($"Количество узлов под 1 краевым после добавления узлов от квадратичного базиса - {FirstConditionNodes.Count}");

            //List<double> Ys = new List<double>();
            
            //for (int i = 0; i < FirstConditionNodes.Count; i++)
            //{
            //    if (Math.Abs(Nodes[FirstConditionNodes[i]].X - _maxX) < 1e-7)
            //        //Console.WriteLine($"X: {Nodes[FirstConditionNodes[i]].X} \t Y: {Nodes[FirstConditionNodes[i]].Y}");
            //        Ys.Add(Nodes[FirstConditionNodes[i]].Y);
            //    if (Math.Abs(Nodes[FirstConditionNodes[i]].X - _maxX) < 1e-7 && Math.Abs(Nodes[FirstConditionNodes[i]].Y - _maxY) < 1e-7)
            //        Console.WriteLine($"X: {Nodes[FirstConditionNodes[i]].X} \t Y: {Nodes[FirstConditionNodes[i]].Y}");
            //}

            //Ys.Sort();
            //Console.WriteLine($"Количество на правой границе : {Ys.Count}");
            //for (int i = 0; i < Ys.Count;i++)
            //{
            //    Console.Write($"{Ys[i]:E5} \t");
            //}
            //Console.WriteLine("\nOUT");
        }
    }

    private void NumberNodes()
    {
        Dictionary<(int, int), (int, int)> hashtable = new Dictionary<(int, int), (int, int)>();
        (int, int) pointInfo;
        int num = Nodes.Count;
        (int, int)[] key = new (int, int)[3];


        // добавляем неосновные узлы для квадратичного базиса
        for (int ielem = 0; ielem < Elements.Length; ielem++)
        {
            key[0] = (Math.Min(Elements[ielem][0], Elements[ielem][1]), Math.Max(Elements[ielem][0], Elements[ielem][1]));
            key[1] = (Math.Min(Elements[ielem][1], Elements[ielem][2]), Math.Max(Elements[ielem][1], Elements[ielem][2]));
            key[2] = (Math.Min(Elements[ielem][0], Elements[ielem][2]), Math.Max(Elements[ielem][0], Elements[ielem][2]));

            for (int i = 0; i < 3; i++)
            {
                if (hashtable.ContainsKey(key[i]))
                {
                    pointInfo = hashtable[key[i]];
                    Elements[ielem][i + 3] = pointInfo.Item1;
                }
                else
                {
                    Elements[ielem][i + 3] = num;
                    pointInfo = (num, ielem);
                    hashtable.Add(key[i], pointInfo);
                    Nodes.Add((Nodes[key[i].Item1] + Nodes[key[i].Item2]) / 2);

                    if (Math.Abs(Nodes[^1].X - _maxX) < 1e-7 ||
                        Math.Abs(Nodes[^1].X - _minX) < 1e-7 ||
                        Math.Abs(Nodes[^1].Y - _maxY) < 1e-7)
                    {
                        FirstConditionNodes.Add(num);
                    }

                    num++;
                }
            }
        }

    }

    private void setMaterials(string muPath, string tokuPath)
    {
        string[][] data1;
        string[][] data2;

        try
        {
            using StreamReader sr1 = new(muPath);
            using StreamReader sr2 = new(tokuPath);

            data1 = sr1.ReadToEnd().Split(new char[] { '\n' }).Select(a => a.Split(' ').Where(str => str != "").ToArray()).ToArray();
            data2 = sr2.ReadToEnd().Split(new char[] { '\n' }).Select(a => a.Split(' ').Where(str => str != "").ToArray()).ToArray();

            if (data1.Length != data2.Length)
                throw new Exception("Количество материалов для токов и мю не одинаково.");

            _materials = new Material[data1.Length - 1];

            for (int i = 0; i < _materials.Length; i++)
            {
                _materials[i] = new Material(Convert.ToDouble(data1[i][1]), Convert.ToDouble(data2[i][1]), data1[i][0]);
            }
        }
        catch (Exception ex)
        {
            Console.WriteLine(ex.Message);
            throw new NotImplementedException();
        }



    }
}
