using FEM2;

Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");
string spaceGridPath = "..\\..\\..\\..\\BinaryConvert\\";
string materialPath = "..\\..\\..\\..\\PRIMER\\";


string splineDataPath = "mu";
(double, double) alphaBeta = (1e-7, 1e-7);
Spline spline = new(200, splineDataPath, alphaBeta);
spline.Compute();
spline.printSpline(10);

bool QuadraticBasis = true;

Grid grid = new(spaceGridPath, materialPath, QuadraticBasis);
FEM fem = new(grid);
fem.SetSpline(spline);
fem.SetNonlinearIterationParametres(20, 1e-12, 1e-14);
fem.SetSlaeParametres(100000, 1e-16);
fem.Compute();
//fem.PrintSolution();

using (StreamWriter sw = new("..\\..\\..\\pointsResults.txt", true))  // append = true
{
    for (int i = 0; i < grid.Materials.Length; i++)
    {
        Material material = grid.Materials[i];
        Console.WriteLine($"{material.Name} {material.Mu} {material.J}");
        sw.WriteLine($"{material.Name} {material.Mu} {material.J:E3}");
    }


    //Console.WriteLine();
    Point2D[] points =
        {
        new Point2D(+5.000000e-002, + 1.100000e-003),
        new Point2D(+ 5.980000e-002, + 3.300000e-003),
        new Point2D(+ 4.000000e-002, + 1.500000e-003),
        new Point2D(+ 4.680000e-002, + 2.100000e-003),
        new Point2D(+ 5.320000e-002, + 2.100000e-003)
    };
    for (int i = 0; i < points.Length; i++)
    {
        (double, double, double, double) res = (fem.AzAtPoint(points[i]), fem.BxAtPoint(points[i]), fem.ByAtPoint(points[i]), fem.AbsBAtPoint(points[i]));
        Console.WriteLine($"Точка ( {points[i].X:E4}; {points[i].Y:E4} ): Az = {res.Item1:.0000E+00}\t  Bx = {res.Item2:.0000E+00}\t  By = {res.Item3:.0000E+00}\t |B| = {res.Item4:.0000E+00}.");
        sw.WriteLine($"Точка ( {points[i].X:E4}; {points[i].Y:E4} ): Az = {res.Item1:.00000000E+00}\t  Bx = {res.Item2:.00000000E+00}\t  By = {res.Item3:.00000000E+00}\t |B| = {res.Item4:.00000000E+00}.");
    }
}

OutputMeshSolution(fem);



void OutputMeshSolution(FEM fem)
{
    double hx = 0.001, hy = 0.001;
    double xStart = -0.25e-2, xEnd = 10.25e-2;
    double yStart = 0.0, yEnd = 6.25e-2;

    int xSteps = Convert.ToInt32((xEnd - xStart) / hx);
    int ySteps = Convert.ToInt32((yEnd - yStart) / hy);

    double x = xStart, y = yStart;

    using StreamWriter sw = new("..\\..\\..\\..\\Graphics\\results.txt");

    for (int i = 0; i <= ySteps; i++)
    {
        x = xStart;

        for (int j = 0; j <= xSteps; j++)
        {
            sw.WriteLine($"{x} {y} {fem.AzAtPoint(new Point2D(x, y))}");
            x = xStart + hx * (j + 1);
        }

        y = yStart + hy * (i + 1);
    }

    Console.WriteLine($"{x}\t {y}");
}