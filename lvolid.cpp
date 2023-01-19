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
    int sock = socket(AF_INET, SOCK_DGRAM, IPPROTO_IP);
    if (sock == -1)
    { /* handle error*/
    }
    ifc.ifc_len = sizeof(buf);
    ifc.ifc_buf = buf;
    if (ioctl(sock, SIOCGIFCONF, &ifc) == -1)
    { /* handle error */
    }
    unsigned char mac_address[6];
    int i = 1;
    for (i = 1; if_indextoname(i, ifr.ifr_name); ++i)
    { // index starts at 1
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
        else
        { /* handle error */
        }
    }
    return 0;
}