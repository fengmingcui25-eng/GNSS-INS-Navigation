#pragma once
#ifndef SOCKETS_H
#define SOCKETS_H

#include <winsock2.h>

bool OpenSocket(SOCKET& sock, const char IP[], const unsigned short Port);
void CloseSocket(SOCKET& sock);

#endif