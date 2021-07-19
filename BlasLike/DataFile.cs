using System;
using System.IO;
namespace DataFile
{
    public class Input
    {
        double[] c = null;
        double[] alpha = null;
        public void Read(string FileName = "testData")
        {
            string line;
            TextReader inF = new StreamReader(FileName);
            while ((line = inF.ReadLine()) != null)
            {
                var lf = line.Split(' ');
                foreach (var k in fields)
                {
                    if (k == "c" && k == lf[0])
                    {
                        c = readDoubleArray(inF);
                    }
                    else if (k == "alpha" && k == lf[0])
                    {
                        alpha = readDoubleArray(inF);
                    }
                }
            }
        }

        string[] fields = "c alpha".Split(' ');
        bool testInFields(int k)
        {
            var back = false;
            if (k == -1) return true;
            foreach (var p in fields)
            {
                if (p[0] == (char)k) back = true;
            }
            return back;
        }
        double[] readDoubleArray(TextReader more)
        {
            string line = more.ReadLine();
            while (!testInFields(more.Peek()))
            {
                line += ' ' + more.ReadLine().Trim();
            }

            double[] back = null;
            var lineF = line.Split(' ');
            var n = 0;
            if ((n = lineF.Length) > 0)
            {
                back = new double[n];
                for (var i = 0; i < n; ++i)
                {
                    back[i] = double.Parse(lineF[i]);
                }
            }
            return back;
        }
    }
}