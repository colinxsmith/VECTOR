#include <sys/ioctl.h>
#include <net/if.h>
#include <unistd.h>
#include <netinet/in.h>
#include <string.h>
#include <cstdio>
int main()
{
    struct ifreq ifr;
    struct ifconf ifc;
    char buf[1028];
    char getname[IF_NAMESIZE];
    int sock = socket(AF_INET, SOCK_DGRAM, IPPROTO_IP);
    if (sock == -1)
    {
        return 1;
    }
    ifc.ifc_len = sizeof(buf);
    ifc.ifc_buf = buf;
    if (ioctl(sock, SIOCGIFCONF, &ifc) == -1)
    {
        return 2;
    }
    unsigned char mac_address[6];
    for (int i = 1; if_indextoname(i, getname); ++i)
    { // index starts at 1
#if defined( __CYGWIN__ )     
    ifr=ifc.ifc_req[i-1];
#else    
       strcpy(ifr.ifr_name, getname);
#endif
        if (ioctl(sock, SIOCGIFFLAGS, &ifr) == 0)
        {
            if (!(ifr.ifr_flags & IFF_LOOPBACK))
            { // don't count loopback
                if (ioctl(sock, SIOCGIFHWADDR, &ifr) == 0)
                {
                    memcpy(mac_address, ifr.ifr_hwaddr.sa_data, 6);
                    printf("Adapter %s:\t%.2x%.2x%.2x%.2x%.2x%.2x\n", ifr.ifr_name, mac_address[0], mac_address[1], mac_address[2], mac_address[3], mac_address[4], mac_address[5]);
                }
            }
        }
    }
    return 0;
}
