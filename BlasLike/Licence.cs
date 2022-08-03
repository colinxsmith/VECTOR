using System;
using System.Text;
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
        public byte[] licenceByteValue = null;
        ///<summary> Write the licence whose data is in licenceByteValue to registry key ourkey </summary>
        ///<param name="ourkey"> string defining registry key </param>
        public bool toRegistry(string ourkey = "Software\\safeqp")
        {
            var back = true;
            if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
            {
                try
                {
                    RegistryKey safekey = Registry.CurrentUser, newkey;
                    if (safekey != null)
                    {
                        newkey = safekey.CreateSubKey(ourkey);
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
                    safekey = Registry.CurrentUser;
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
                    Console.WriteLine($"Newly Written Licence: \t{licence}");
                }
                catch (Exception prob)
                {
                    back = false;
                    Console.WriteLine("exception" + prob);
                }
            }
            return back;
        }
        ///<summary> Read the licence in registry key ourkey to licenceByteValue  </summary>
        ///<param name="ourkey"> string defining registry key </param>
        public bool fromRegistry(string ourkey = "Software\\safeqp")
        {
            bool worked = true;
            if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
            {
                worked = false;
                try
                {
                    RegistryKey safekey = Registry.CurrentUser, newkey;
                    if (safekey != null)
                    {
                        newkey = safekey.OpenSubKey(ourkey);
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
                    ColourConsole.WriteEmbeddedColourLine($"[green]Our Key =[/green] \t\t[yellow]{ourkey}[/yellow]");
                    ColourConsole.WriteEmbeddedColourLine($"[green]Current Licence:[/green] \t[yellow]{licence}[/yellow]");
                    safekey.Dispose();
                }
                catch (Exception prob)
                {
                    ColourConsole.WriteError("exception" + prob);
                }
            }
            return worked;
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
            System.DateTime dtDateTime = new DateTime(1970, 1, 1, 0, 0, 0, 0, System.DateTimeKind.Utc);
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
        public UInt32 VolId()
        {
            var hex = "facc0ff5";//Start with this in case there are no more
            var output = Convert.ToUInt32(hex, 16);
            ColourConsole.WriteEmbeddedColourLine($"[green]{hex} =[/green] [cyan]{output}[/cyan] [yellow]{output:x}[/yellow]");
            if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
                foreach (var nic in NetworkInterface.GetAllNetworkInterfaces())
                {
                    if (nic.NetworkInterfaceType.ToString().Contains("Ethernet"))//&& nic.NetworkInterfaceType.ToString().Contains("USB"))
                    {
                        if (nic.Description.ToString().Contains("Virtual")) continue;
                        if (nic.Description.ToString().Contains("USB")) continue;
                        /*   Console.WriteLine(nic.Name);
                           Console.WriteLine(nic.GetType());
                           Console.WriteLine(nic.NetworkInterfaceType);
                           Console.WriteLine(nic.GetPhysicalAddress());
                           Console.WriteLine();*/
                        var address = nic.GetPhysicalAddress().ToString();
                        ColourConsole.WriteInfo($"Adapter: {address}");
                        var bAddress = new UInt32[6];
                        for (var i = 0; i < 6; i++)
                        {
                            bAddress[i] = Convert.ToUInt32(address.Substring(i * 2, 2), 16);
                            ColourConsole.Write($"{bAddress[i]} ");
                        }
                        Array.Sort(bAddress);
                        Array.Reverse(bAddress);
                        for (var i = 0; i < 6; i++) ColourConsole.Write($"{bAddress[i],2:x} ", ConsoleColor.DarkMagenta);
                        Console.WriteLine();
                        var newByte = new UInt32[4];
                        newByte[0] = bAddress[0];
                        newByte[1] = bAddress[1];
                        newByte[2] = ((bAddress[2] + bAddress[3]) / 2) & 0xff;
                        newByte[3] = ((bAddress[4] + bAddress[5]) / 2) & 0xff;
                        var bStr = $"{newByte[0],2:x2}{newByte[1],2:x2}{newByte[2],2:x2}{newByte[3],2:x2}";
                        ColourConsole.WriteInfo(bStr);
                        bStr = bStr.Replace(" ", "");
                        var bStrInt = Convert.ToUInt32(bStr, 16);
                        ColourConsole.WriteEmbeddedColourLine($"[green]{bStr} =[/green] [cyan]{bStrInt}[/cyan] [yellow]{bStrInt:x}[/yellow]");
                        output ^= bStrInt;
                    }
                }
            ColourConsole.WriteEmbeddedColourLine($"Finally [cyan]{output}[/cyan] [yellow]{output:x}[/yellow]");
            return output;
        }
        public string VersionString(){
            const string back= "Optimiser from BITA Plus. \0aaaabbbbccccddddeeeeffffgggghhhhiiiiJJJJ";
            return back;
        }
    }
}
