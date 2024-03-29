using System;
using System.IO;
using System.Collections.Generic;
namespace DataFile
{
    public class InputSomeData : IDisposable
    {
        public void Dispose()
        {
            if (mapDouble != null)
            {
                foreach (var dv in mapDouble.Keys)
                {
                    mapDouble[dv] = null;
                }
                mapDouble.Clear();
            }
            if (mapInt != null)
            {
                foreach (var dv in mapInt.Keys)
                {
                    mapInt[dv] = null;
                }
                mapInt.Clear();
            }
            if (mapString != null)
            {
                foreach (var dv in mapString.Keys)
                {
                    mapString[dv] = null;
                }
                mapString.Clear();
            }
        }
        public Dictionary<string, double[]> mapDouble = null;
        public Dictionary<string, int[]> mapInt = null;
        public Dictionary<string, string[]> mapString = null;
        public char dataSep = ' ';
        public string doubleFields = "c alpha";
        public string intFields = "n m";
        public string stringFields = "names";
        public void Write(string OutFile = "out.log", bool append = false)
        {
            using (var outF = new StreamWriter(OutFile, append))
            {
                foreach (var name in intFields.Split(dataSep))
                {
                    outF.WriteLine(name);
                    var n = mapInt[name].Length - 1;
                    for (var i = 0; i < n; ++i)
                    {
                        outF.Write(mapInt[name][i] + " ");
                    }
                    outF.WriteLine(mapInt[name][n]);
                }
                foreach (var name in doubleFields.Split(dataSep))
                {
                    outF.WriteLine(name);
                    var n = mapDouble[name].Length - 1;
                    for (var i = 0; i < n; ++i)
                    {
                        outF.Write(mapDouble[name][i] + " ");
                    }
                    outF.WriteLine(mapDouble[name][n]);
                }
                foreach (var name in stringFields.Split(dataSep))
                {
                    outF.WriteLine(name);
                    var n = mapString[name].Length - 1;
                    for (var i = 0; i < n; ++i)
                    {
                        outF.Write(mapString[name][i] + " ");
                    }
                    outF.WriteLine(mapString[name][n]);
                }
            }
        }
        public void Read(string FileName = "testData")
        {
            using (var inF = new StreamReader(FileName))
            {
                string line;
                bool setData = false;
                while ((line = inF.ReadLine()) != null)
                {if(line.Contains("---------"))break;
                    line = line.Trim();
                    var lf = line.Split(dataSep);
                    setData = false;
                    if (line == "") continue;
                    if (!setData)
                        foreach (var k in doubleFields.Split(dataSep))
                        {
                            if (k == lf[0])
                            {
                                if (mapDouble == null) mapDouble = new Dictionary<string, double[]>();
                                mapDouble[k] = readDoubleArray(inF);
                                setData = true; break;
                            }
                        }
                    if (!setData)
                        foreach (var k in intFields.Split(dataSep))
                        {
                            if (k == lf[0])
                            {
                                if (mapInt == null) mapInt = new Dictionary<string, int[]>();
                                mapInt[k] = readIntArray(inF);
                                setData = true; break;
                            }
                        }
                    if (!setData)
                        foreach (var k in stringFields.Split(dataSep))
                        {
                            if (k == lf[0])
                            {
                                if (mapString == null) mapString = new Dictionary<string, string[]>();
                                mapString[k] = readStringArray(inF);
                                setData = true; break;
                            }
                        }
                }
            }
        }

        bool testInFields(int k)
        {
            var back = false;
            if (k == -1) return true;
            if (doubleFields != "")
                foreach (var p in doubleFields.Split(dataSep))
                {
                    if (p[0] == (char)k) back = true;
                }
            if (intFields != "")
                foreach (var p in intFields.Split(dataSep))
                {
                    if (p[0] == (char)k) back = true;
                }
            if (stringFields != "")
                foreach (var p in stringFields.Split(dataSep))
                {
                    if (p[0] == (char)k) back = true;
                }
            return back;
        }
        // Check that a line of entries "out of place" contains doubles or integers.
        // so that it can be added on to the previous line. We assume all doubles, integers or strings on a line.
        // This should handle most mistakes in a data file.
        double[] readDoubleArray(StreamReader more)
        {
            var line = more.ReadLine();
            bool foundDouble = true;
            if (line != null) line = line.Trim();
            while (!testInFields(more.Peek()))
            {
                var ll = more.ReadLine();
                try
                {
                    if (ll != null && !ll.Contains(dataSep)) double.Parse(ll);
                    else if (ll != null && ll.Contains(dataSep)) double.Parse(ll.Split(dataSep)[0]);
                }
                catch
                {
                    foundDouble = false;
                }
                if (foundDouble && ll != null && ll != "") line += ' ' + ll.Trim();
            }
            double[] back = null;
            var lineF = line.Split(dataSep);
            var n = 0;
            if (line != "" && (n = lineF.Length) > 0)
            {
                back = new double[n];
                for (var i = 0; i < n; ++i)
                {
                    back[i] = double.Parse(lineF[i]);
                }
            }
            return back;
        }
        int[] readIntArray(StreamReader more)
        {
            var line = more.ReadLine();
            bool foundInt = true;
            if (line != null) line = line.Trim();
            while (!testInFields(more.Peek()))
            {
                var ll = more.ReadLine(); try
                {
                    if (ll != null && !ll.Contains(dataSep)) int.Parse(ll);
                    if (ll != null && ll.Contains(dataSep)) int.Parse(ll.Split(dataSep)[0]);
                }
                catch
                {
                    foundInt = false;
                }
                if (foundInt && ll != null) line += ' ' + ll.Trim();
            }
            int[] back = null;
            var lineF = line.Split(dataSep);
            var n = 0;
            if (line != "" && (n = lineF.Length) > 0)
            {
                back = new int[n];
                for (var i = 0; i < n; ++i)
                {
                    back[i] = int.Parse(lineF[i]);
                }
            }
            return back;
        }
        string[] readStringArray(StreamReader more)
        {
            var line = more.ReadLine();
            if (line != null) line = line.Trim();
            while (!testInFields(more.Peek()))
            {
                var ll = more.ReadLine();
                if (ll != null) line += ' ' + ll.Trim();
            }
            string[] back = null;
            var lineF = line.Split(dataSep);
            var n = 0;
            if ((n = lineF.Length) > 0)
            {
                back = new string[n];
                for (var i = 0; i < n; ++i)
                {
                    back[i] = lineF[i];
                }
            }
            return back;
        }
        public void PrintField(string n)
        {
            if (mapString.ContainsKey(n))
                ActiveSet.Optimise.printV(n, mapString[n]);
            else if (mapInt.ContainsKey(n))
                ActiveSet.Optimise.printV(n, mapInt[n]);
            else if (mapDouble.ContainsKey(n))
                ActiveSet.Optimise.printV(n, mapDouble[n]);
        }
    }
}