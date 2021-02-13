/* fpont 12/99 */
/* pont.net    */
/* udpserver.c */

/* Converted to echo client/server with select() (timeout option).
   See testTUDPSocket.C */
/* Compile with: gcc udpserver.c -o udpserver */
/* on Windows: cl -nologo -Z7 -MD -GR -EHsc udpserver.c */
/* 3/30/05 John Schultz */

#include <stdlib.h>
#include <sys/types.h>
#ifdef _WIN32
#include <winsock2.h>
#include <windows.h>
#else
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <unistd.h> /* close() */
#endif
#include <stdio.h>
#include <string.h> /* memset() */

#ifdef _WIN32
#pragma comment(lib,"Ws2_32.lib")
#endif

#define LOCAL_SERVER_PORT 1500
#define MAX_MSG 100

int main(int argc, char *argv[]) {

   int sd, rc, n, flags;
   unsigned cliLen;
   struct sockaddr_in cliAddr, servAddr;
   char msg[MAX_MSG];

#ifdef _WIN32
   WSADATA wsaData = {0};
   int iResult = 0;
   /* Initialize Winsock */
   iResult = WSAStartup(MAKEWORD(2, 2), &wsaData);
   if (iResult != 0) {
      wprintf(L"WSAStartup failed: %d\n", iResult);
      return -1;
   }
   /* socket creation */
   sd=socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
#else
   /* socket creation */
   sd=socket(AF_INET, SOCK_DGRAM, 0);
#endif
   if (sd<0) {
      printf("%s: cannot open socket \n",argv[0]);
      exit(1);
   }

   /* bind local server port */
   servAddr.sin_family = AF_INET;
   servAddr.sin_addr.s_addr = htonl(INADDR_ANY);
   servAddr.sin_port = htons(LOCAL_SERVER_PORT);
   rc = bind (sd, (struct sockaddr *) &servAddr, sizeof(servAddr));
   if (rc<0) {
      printf("%s: cannot bind port number %d \n", argv[0], LOCAL_SERVER_PORT);
      exit(1);
   }
   printf("%s: waiting for data on port UDP %u\n", argv[0], LOCAL_SERVER_PORT);

/* BEGIN jcs 3/30/05 */
   flags = 0;
/* END jcs 3/30/05 */

   /* server infinite loop */
   while (1) {

      /* init buffer */
      memset(msg,0x0,MAX_MSG);

      /* receive message */
      cliLen = sizeof(cliAddr);
      n = recvfrom(sd, msg, MAX_MSG, flags, (struct sockaddr *) &cliAddr,
                   &cliLen);

      if (n<0) {
         printf("%s: cannot receive data \n",argv[0]);
         continue;
      }
      /* print received message */
      printf("%s: from %s:UDP%u : %s \n", argv[0],
             inet_ntoa(cliAddr.sin_addr), ntohs(cliAddr.sin_port), msg);

/* BEGIN jcs 3/30/05 */

#ifdef _WIN32
      Sleep(1000);
#else
      sleep(1);
#endif
      sendto(sd, msg, n, flags, (struct sockaddr *)&cliAddr, cliLen);

/* END jcs 3/30/05 */

   }/* end of server infinite loop */
   return 0;
}
