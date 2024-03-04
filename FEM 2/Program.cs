using UMFCourseProject;

Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");
string spaceGridPath = "C:\\Users\\Skromer\\source\\repos\\FEM 2\\BinaryConvert\\";
string materialPath = "C:\\TELMA\\PRIMER\\";
Grid grid = new(spaceGridPath, materialPath);
FEM fem = new(grid);
fem.Compute();
//fem.PrintSolution();
for (int i = 0; i < grid.Materials.Length; i++)
{
    Material material = grid.Materials[i];
    Console.WriteLine($"{material.Name} {material.Mu} {material.J}");
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
    double res = fem.ValueAtPoint(points[i]);
    Console.WriteLine($"Значение вектор-потенциала в точке ({points[i].X};{points[i].Y}) равно {res:E11}.");
}
