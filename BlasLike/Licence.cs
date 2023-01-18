using System;
using System.Text;
using System.IO;
using System.Net.NetworkInformation;
using Microsoft.Win32;
using System.Runtime.InteropServices;
using Microsoft.Extensions.Hosting.WindowsServices;
namespace Licensing
{
    [StructLayout(LayoutKind.Explicit)]
    public struct byteint
    {//1 int has length 4 bytes
        [FieldOffset(0)]
        public int mainint;
        [FieldOffset(0)]
        public byte byte1;
        [FieldOffset(1)]
        public byte byte2;
        [FieldOffset(2)]
        public byte byte3;
        [FieldOffset(3)]
        public byte byte4;
    }
    public class validator_t
    {
        public validator_t Copy()
        {
            var back = new validator_t();
            back.b = (byte[])this.b.Clone();
            return back;
        }
        public byte[] b = new byte[16];
        public int pad
        {
            get
            {
                byteint conv = new byteint();
                conv.mainint = 0;
                conv.byte1 = b[0];
                conv.byte2 = b[1];
                conv.byte3 = b[2];
                conv.byte4 = b[3];
                return conv.mainint;
            }
            set
            {
                byteint conv = new byteint();
                conv.mainint = value;
                b[0] = conv.byte1;
                b[1] = conv.byte2;
                b[2] = conv.byte3;
                b[3] = conv.byte4;
            }
        }
        public int start
        {
            get
            {
                byteint conv = new byteint();
                conv.mainint = 0;
                conv.byte1 = b[4];
                conv.byte2 = b[5];
                conv.byte3 = b[6];
                conv.byte4 = b[7];
                return conv.mainint;
            }
            set
            {
                byteint conv = new byteint();
                conv.mainint = value;
                b[4] = conv.byte1;
                b[5] = conv.byte2;
                b[6] = conv.byte3;
                b[7] = conv.byte4;
            }
        }
        public int stop
        {
            get
            {
                byteint conv = new byteint();
                conv.mainint = 0;
                conv.byte1 = b[8];
                conv.byte2 = b[9];
                conv.byte3 = b[10];
                conv.byte4 = b[11];
                return conv.mainint;
            }
            set
            {
                byteint conv = new byteint();
                conv.mainint = value;
                b[8] = conv.byte1;
                b[9] = conv.byte2;
                b[10] = conv.byte3;
                b[11] = conv.byte4;
            }
        }
        public int hid
        {
            get
            {
                byteint conv = new byteint();
                conv.mainint = 0;
                conv.byte1 = b[12];
                conv.byte2 = b[13];
                conv.byte3 = b[14];
                conv.byte4 = b[15];
                return conv.mainint;
            }
            set
            {
                byteint conv = new byteint();
                conv.mainint = value;
                b[12] = conv.byte1;
                b[13] = conv.byte2;
                b[14] = conv.byte3;
                b[15] = conv.byte4;
            }
        }
    }
    public class Licence
    {
        int bitaopt = 23;
        public string VersionString = "";
        public Licence()
        {
            validator_m_byte = new byte[validator_m.Length];
            for (var i = 0; i < validator_m.Length; ++i)
            {
                validator_m_byte[i] = (byte)validator_m[i];
            }
        }
        public string validator_m = "RdG^AsTi1GgaKa7teDF3E645Fd912j4cB4b0_1|(&50'34$!TzZjkIXlMHPVW*@~";
        public byte[] validator_m_byte;
        byte validator_c(int z, int k) => validator_m_byte[(z + k) % 48];
        string licence = "";
       public string connectedNames="";
        public byte[] licenceByteValue = null;
        public bool deleteKey(string ourkey = "Software\\safeqp")
        {
            RegistryKey safekey, newkey;
            if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
            {
                try
                {
                    safekey = Registry.LocalMachine; newkey = safekey.OpenSubKey(ourkey, true);
                    if (newkey == null)
                    {
                        newkey = safekey.CreateSubKey(ourkey, true);
                        newkey.DeleteSubKey(ourkey);
                    }
                    safekey.Dispose();
                }
                catch { safekey = Registry.CurrentUser; }
                newkey = safekey.CreateSubKey(ourkey, true);
                try
                {
                    newkey.DeleteValue(ourkey);
                    safekey.Dispose();
                }
                catch { return false; }
            }
            else
            {
                var basef = AppContext.BaseDirectory + "licence";
                File.Delete(basef);
            }
            return true;
        }
        ///<summary> Write the licence whose data is in licenceByteValue to registry key ourkey </summary>
        ///<param name="ourkey"> string defining registry key </param>
        public bool toRegistry(bool usefile = false, string ourkey = "Software\\safeqp", bool print = false)
        {
            var back = true;
            if (!usefile && RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
            {
                try
                {
                    RegistryKey safekey, newkey;
                    try
                    {
                        safekey = Registry.LocalMachine; newkey = safekey.OpenSubKey(ourkey, true);
                        safekey.Dispose();
                    }
                    catch { safekey = Registry.CurrentUser; }
                    if (safekey != null)
                    {
                        newkey = safekey.CreateSubKey(ourkey, true);
                        if (newkey == null)
                        {
                            safekey.Dispose(); back = false;
                            throw new Exception("No key");
                        }
                        if (licenceByteValue != null && licenceByteValue.Length > 0)
                            newkey.SetValue(ourkey, licenceByteValue);
                        else if (newkey.GetValue(ourkey) != null)
                            newkey.DeleteValue(ourkey);
                    }
                    safekey.Dispose();
                    try
                    {
                        safekey = Registry.LocalMachine; newkey = safekey.OpenSubKey(ourkey, true);
                        safekey.Dispose();
                    }
                    catch { safekey = Registry.CurrentUser; }
                    newkey = safekey.OpenSubKey(ourkey);
                    licenceByteValue = (Byte[])newkey.GetValue(ourkey);
                    if (licenceByteValue == null)
                    {
                        safekey.Dispose(); back = false;
                        throw new Exception("No licence!");
                    }
                    safekey.Dispose();
                    licence = "";
                    for (int i = 0; i < licenceByteValue.Length; ++i)
                    {
                        licence += string.Format("{0:x2};", licenceByteValue[i]);
                    }
                    if (print) Console.WriteLine($"Newly Written Licence: \t{licence}");
                }
                catch (Exception prob)
                {
                    back = false;
                    Console.WriteLine("exception" + prob);
                }
            }
            else
            {
                var basef = AppContext.BaseDirectory + "licence";
                if (print) ColourConsole.WriteInfo(basef);
                try
                {
                    using (var stream = File.Open(basef, FileMode.OpenOrCreate))
                    {
                        using (var write = new BinaryWriter(stream, Encoding.UTF8))
                        {
                            for (var i = 0; i < licenceByteValue.Length; ++i)
                                write.Write(licenceByteValue[i]);
                        }
                    }
                }
                catch { back = false; }
            }
            return back;
        }
        ///<summary> Read the licence in registry key ourkey to licenceByteValue 
        ///Returns; 0 if failed, 1 if run as root, 2 if run as user </summary>
        ///<param name="ourkey"> string defining registry key </param>
        public int fromRegistry(bool usefile = false, string ourkey = "Software\\safeqp", bool print = false)
        {
            bool worked = true;
            bool root = true;
            if (!usefile && RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
            {
                worked = false;
              //  try
                {
                    RegistryKey safekey, newkey;
                    try
                    {
                        safekey = Registry.LocalMachine; newkey = safekey.OpenSubKey(ourkey, true);
                        safekey.Dispose();
                    }
                    catch { safekey = Registry.CurrentUser; root = false; }
                    if (safekey != null)
                    {
                        newkey = safekey.OpenSubKey(ourkey, true);
                        if (newkey == null)
                        {
                            safekey.Dispose();
                            throw new Exception("No key");
                        }
                        licenceByteValue = (Byte[])newkey.GetValue(ourkey);
                        if (licenceByteValue == null)
                        {
                            safekey.Dispose();
                            throw new Exception("No licenceByteValue");
                        }
                        else worked = true;
                        licence = "";
                        for (int i = 0; i < licenceByteValue.Length; ++i)
                        {
                            licence += $"{licenceByteValue[i],2:x2};";
                        }
                    }
                    if (print) ColourConsole.WriteEmbeddedColourLine($"[green]Our Key =[/green] \t\t[yellow]{ourkey}[/yellow]");
                    if (print) ColourConsole.WriteEmbeddedColourLine($"[green]Current Licence:[/green] \t[yellow]{licence}[/yellow]");
                    safekey.Dispose();
                }
         /*       catch (Exception prob)
                {
                    ColourConsole.WriteError("exception" + prob);
                }*/
            }
            else
            {
                worked = false;
                var basef = AppContext.BaseDirectory + "licence";
                if (print) ColourConsole.WriteInfo(basef);
                try
                {
                    using (var stream = File.Open(basef, FileMode.Open))
                    {
                        using (var read = new BinaryReader(stream, Encoding.UTF8, false))
                        {
                            licenceByteValue = read.ReadBytes(20);
                        }
                        worked = true; root = false;
                    }
                }
                catch {; }
            }
            if (!worked) return 0;
            return worked && root ? 1 : 2;
        }
        public void krypton(byte[] b = null, int bstart = 0)
        {
            if (b == null) b = licenceByteValue;
            int i;
            byte[] w = new byte[16], z = new byte[16];
            for (i = 0; i < 16; i++)
            {
                w[i] = (byte)(b[i + bstart] >> (byte)4);
                z[i] = (byte)(b[i + bstart] & (byte)0xf);
            }
            for (i = 0; i < 16; i++)
            {
                b[i + bstart] = (byte)((z[15 - i] << (byte)4) | (w[15 - i]));
            }
        }
        public void make_valid(validator_t vp, int start, int stop, int hid)
        {
            int i, j, k;
            byte c;

            for (j = 0; j < 48; j += 16) krypton(validator_m_byte, j);
            vp.pad = 0x13101955 + bitaopt;
            vp.start = start;
            vp.stop = stop;
            vp.hid = hid;
            for (i = j = 0; i < 16; i++) j += vp.b[i];
            k = (j * 147) % 48;

            vp.pad += stop - start;
            vp.pad ^= hid;
            vp.start ^= vp.pad;
            vp.hid ^= stop;

            krypton(vp.b);
            for (j = 0; j < 3; j++)
            {
                for (i = 0; i < 16; i++) vp.b[i] ^= validator_c(i + j * 16, k);
                krypton(vp.b);
                for (c = vp.b[i = 0]; i < 15; i++) vp.b[i] = vp.b[i + 1];
                vp.b[i] = c;
            }
            for (i = 0; i < 16; i++) vp.b[i] ^= validator_c(i + j * 16, k);

            /*	printf((char*)"%8.8lX %8.8lX %8.8lX %8.8lX\n", */
            /*		vp->t.start,vp->t.stop,vp->t.pad,vp->t.hid); */
            for (j = 0; j < 48; j += 16) krypton(validator_m_byte, j);
        }
        public DateTime UnixTimeStampToDateTime(double unixTimeStamp)
        {
            // Unix timestamp is seconds past epoch
            DateTime dtDateTime = new DateTime(1970, 1, 1, 0, 0, 0, 0, System.DateTimeKind.Utc);
            dtDateTime = dtDateTime.AddSeconds(unixTimeStamp).ToLocalTime();
            return dtDateTime;
        }
        public bool check_valid(ref validator_t vp)
        {
            int i, j, k;
            byte c;
            validator_t v = vp.Copy();

            for (j = 0; j < 48; j += 16) krypton(validator_m_byte, j);
            for (k = 0; k < 48; k++)
            {
                v = vp.Copy();
                //	printf((char*)"k=%d\n",k);
                //	printf((char*)"%8.8lX %8.8lX %8.8lX %8.8lX\n",
                //		v.t.start,v.t.stop,v.t.pad,v.t.hid);
                for (j = 3, i = 0; i < 16; i++) v.b[i] ^= validator_c(i + j * 16, k);
                while (j-- > 0)
                {
                    for (c = v.b[i = 15]; i > 0; i--) v.b[i] = v.b[i - 1];
                    v.b[i] = c;
                    krypton(v.b);
                    for (i = 0; i < 16; i++) v.b[i] ^= validator_c(i + j * 16, k);
                }

                krypton(v.b);
                //	printf((char*)"%8.8lX %8.8lX %8.8lX %8.8lX\n",
                //		v.t.start,v.t.stop,v.t.pad,v.t.hid);
                v.hid ^= v.stop;
                v.start ^= v.pad;
                v.pad ^= v.hid;
                v.pad -= v.stop - v.start;

                //	printf((char*)"%8.8lX %8.8lX %8.8lX %8.8lX\n",
                //		v.t.start,v.t.stop,v.t.pad,v.t.hid);
                if (v.pad == 0x13101955 + bitaopt)
                {
                    vp = v;
                    break;
                }
            }
            for (j = 0; j < 48; j += 16) krypton(validator_m_byte, j);
            //	std::cout<< "compare "<<strcmp(validator_m,vv)<<std::endl;

            return k != 48;
        }
        ///<summary>If hid is 0, decode data in s to find hid,start and stop,
        ///if hid is not 0, encode hid, start and stop to find s.</summary>
        ///<param name="s">byte array (length 16) of encoded licence data</param>
        ///<param name="hid">Volume id of computer (or 0x13101955 for set up)</param>
        ///<param name="start">integer defining the licence start (current time is now.ToUnixTimeSeconds())</param>
        ///<param name="stop">integer defining the licence stop (current time is now.ToUnixTimeSeconds())</param>
        public void convert(byte[] s, ref int hid, ref int start, ref int stop)
        {//Helper for make_valid and check_valid
            validator_t validator = new validator_t();
            if (hid != 0)
            {
                make_valid(validator, start, stop, hid);
                for (var i = 0; i < 16; ++i) s[i] = validator.b[i];
            }
            else
            {
                for (var i = 0; i < 16; ++i) validator.b[i] = s[i];
                check_valid(ref validator);
                hid = validator.hid;
                start = validator.start;
                stop = validator.stop;
            }
        }
        public Int32 VolId(bool usefile = false, bool print = false)
        {
            var hex = "facc0ff5";//Start with this in case there are no more
            var output = Convert.ToInt32(hex, 16);
            if (print) ColourConsole.WriteEmbeddedColourLine($"[green]{hex} =[/green] [cyan]{output}[/cyan] [yellow]{output:x}[/yellow]");
            //if (!usefile&&RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
            foreach (var nic in NetworkInterface.GetAllNetworkInterfaces()) //This works on linux
            {
                if (nic.NetworkInterfaceType.ToString().Contains("Ethernet")||nic.NetworkInterfaceType.ToString().Contains("Wireless"))//&& nic.NetworkInterfaceType.ToString().Contains("USB"))
                {
                    if (nic.Description.ToString().ToLower().Contains("virtual")) continue;
                    if (nic.Description.ToString().ToLower().Contains("usb")) continue;
                    if (nic.Description.ToString().ToLower().Contains("bluetooth")) continue;
                    /*   Console.WriteLine(nic.Name);
                       Console.WriteLine(nic.GetType());
                       Console.WriteLine(nic.NetworkInterfaceType);
                       Console.WriteLine(nic.GetPhysicalAddress());
                       Console.WriteLine();*/
                       connectedNames+=nic.Description+";";
                    var address = nic.GetPhysicalAddress().ToString();
                    if (print) ColourConsole.WriteInfo($"Adapter: {address}");
                    var bAddress = new int[6];
                    for (var i = 0; i < 6; i++)
                    {
                        bAddress[i] = Convert.ToInt16(address.Substring(i * 2, 2), 16);
                        if (print) ColourConsole.Write($"{bAddress[i]} ");
                    }
                    Array.Sort(bAddress);
                    Array.Reverse(bAddress);
                    if (print)
                    {
                        for (var i = 0; i < 6; i++) ColourConsole.Write($"{bAddress[i],2:x} ", ConsoleColor.DarkMagenta);
                        Console.WriteLine();
                    }
                    var newByte = new Int32[4];
                    newByte[0] = bAddress[0];
                    newByte[1] = bAddress[1];
                    newByte[2] = ((bAddress[2] + bAddress[3]) / 2) & 0xff;
                    newByte[3] = ((bAddress[4] + bAddress[5]) / 2) & 0xff;
                    var bStr = $"{newByte[0],2:x2}{newByte[1],2:x2}{newByte[2],2:x2}{newByte[3],2:x2}";
                    if (print) ColourConsole.WriteInfo(bStr);
                    bStr = bStr.Replace(" ", "");
                    var bStrInt = Convert.ToInt32(bStr, 16);
                    if (print) ColourConsole.WriteEmbeddedColourLine($"[green]{bStr} =[/green] [cyan]{bStrInt}[/cyan] [yellow]{bStrInt:x}[/yellow]");
                    output ^= bStrInt;
                }
            }
            if (print) ColourConsole.WriteEmbeddedColourLine($"Finally [cyan]{output}[/cyan] [yellow]{output:x}[/yellow]");
            return output;
        }
        public bool CheckLicence(bool print = false, bool usefile = false)
        {
            var rootPath = AppContext.BaseDirectory;
            DateTime fileTime = File.GetLastWriteTime(rootPath + "BlasLike.dll");
            const string version = "1.0";
            var back = $"BITA Plus ASP.NET Core Portfolio Optimiser Version {version} made {fileTime}";
            var pass = false;
            var vid = VolId(usefile);
            int start = 0, stop = 0, hid = 0;
            string printStart = "", printStop = "", printNow = "";
            int fromReg;
            try{
            if ((fromReg = fromRegistry(usefile)) > 0)
            {
                DateTimeOffset now = new DateTimeOffset(DateTime.Now);
                var year = now.Year;
                var month = now.Month;
                var day = now.Day;
                var now1 = new DateTime(year, month, day, 1, 2, 3, 4);
                var now11 = new DateTimeOffset(now1);
                var newstart = now11.ToUnixTimeSeconds();
                printNow = $"{now}";
                var timenow = now.ToUnixTimeSeconds();
                if (timenow < newstart) newstart -= 24 + 2 * 60 + 3;
                convert(licenceByteValue, ref hid, ref start, ref stop);
                printStart = $"{UnixTimeStampToDateTime(start)}";
                printStop = $"{UnixTimeStampToDateTime(stop)}";
                Licensing.byteint curveKeys = new Licensing.byteint();
                curveKeys.byte1 = licenceByteValue[16];
                curveKeys.byte2 = licenceByteValue[17];
                curveKeys.byte3 = licenceByteValue[18];
                curveKeys.byte4 = licenceByteValue[19];
                var ckeys = Convert.ToString(curveKeys.mainint, 2);
                hid -= curveKeys.mainint;
                pass = true;
                pass = pass && (start < timenow);
                pass = pass && (stop > timenow);
                pass = pass && (hid == vid || hid == 0x13101955);
                var badhid=!(hid == vid || hid == 0x13101955);
                var extramess=badhid?$"hid should be {vid:x}":"";
                string user = fromReg == 1 ? "root" : "user";
                int days = (stop - (int)timenow) / 24 / 3600;
                if (pass) back += $".\nRunning as {user}. Licence starts: {printStart} until: {printStop}. I.e. {days} days left.\nTime now: {printNow}. Valid on: {hid:x}.\nKeys: {ckeys}";
                else back += $".\nRunning as {user}. Licence is not valid!!!!!!!!!!!!! From: {printStart} until: {printStop}.\nTime now: {printNow}. Valid on: {hid:x}.\nKeys: {ckeys} {extramess}";
                if (pass)
                {//Reset the start time and change hid to the that for this machine
                    start = (int)newstart;
                    hid = (int)vid;
                    hid += curveKeys.mainint;
                    convert(licenceByteValue, ref hid, ref start, ref stop);
                    toRegistry(usefile);
                }
                else deleteKey();
            }
            else{
                back+=$" No licence found!!!!";
                pass=false;
            }}
            catch{
                back+=$" No licence found!!!!";
                pass=false;
            }
            if (print) ColourConsole.WriteEmbeddedColourLine(back);
            VersionString = back;
            return pass;
        }
    }
}
