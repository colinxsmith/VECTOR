using System;
using Microsoft.Win32;
using System.Runtime.InteropServices;
namespace Licensing
{
    [StructLayout(LayoutKind.Explicit)]
    struct byteint
    {
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
        public byte[] b;
        public int pad
        {byteint conv = new byteint(); 
            get
            {
                conv.mainint = 0;
                conv.byte1 = b[0];
                conv.byte2 = b[1];
                conv.byte3 = b[2];
                conv.byte4 = b[3];
                return conv.mainint;
            }
            set
            {
                conv.mainint = value;
                b[0] = conv.byte1;
                b[1] = conv.byte2;
                b[2] = conv.byte3;
                b[3] = conv.byte4;
            }
        }
        public int start
        {byteint conv = new byteint();
            get
            {
                conv.mainint = 0;
                conv.byte1 = b[4];
                conv.byte2 = b[5];
                conv.byte3 = b[6];
                conv.byte4 = b[7];
                return conv.mainint;
            }
            set
            {
                conv.mainint = value;
                b[4] = conv.byte1;
                b[5] = conv.byte2;
                b[6] = conv.byte3;
                b[7] = conv.byte4;
            }
        }
        public int stop
        {byteint conv = new byteint();
            get
            {
                conv.mainint = 0;
                conv.byte1 = b[8];
                conv.byte2 = b[9];
                conv.byte3 = b[10];
                conv.byte4 = b[11];
                return conv.mainint;
            }
            set
            {
                conv.mainint = value;
                b[8] = conv.byte1;
                b[9] = conv.byte2;
                b[10] = conv.byte3;
                b[11] = conv.byte4;
            }
        }
        public int hid
        {byteint conv = new byteint();
            get
            {
                conv.mainint = 0;
                conv.byte1 = b[12];
                conv.byte2 = b[13];
                conv.byte3 = b[14];
                conv.byte4 = b[15];
                return conv.mainint;
            }
            set
            {
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
        public byte[] licenceByteValue;
        public bool fromRegistry(string ourkey = "Software\\safeqp")
        {
            bool worked = false;
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
        public bool check_valid(ref validator_t vp)
        {
            int i, j, k;
            byte c;
            validator_t v;

            for (j = 0; j < 48; j += 16) krypton(validator_m_byte, j);
            for (k = 0; k < 48; k++)
            {
                v = vp;
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
    }
}
