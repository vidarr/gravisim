/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*- */
/*
 * viewer.c
 * Copyright (C) Michael J. Beer 2014 <michael@ubeer.org>
 * 
 * uniview is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * uniview is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/*--------------------------------------------------------------------------*/
#include <SDL/SDL.h>
#include <stdio.h>
#include <time.h>
/*--------------------------------------------------------------------------*/
static int xResolution = 640;
static int yResolution = 480;
static int framesPerSecond = 1;
/*--------------------------------------------------------------------------*/
#define GOT_NEW_VECTOR 1
/*----------------------------------------------------------------------------*/
void setPixel(SDL_Surface *surface, int x, int y, Uint32 pixel)
    {
        int bpp = surface->format->BytesPerPixel;
        /* Here p is the address to the pixel we want to set */
        Uint8 *p = (Uint8 *)surface->pixels + y * surface->pitch + x * bpp;
        switch(bpp) {
            case 1:
                *p = pixel;
                break;
            case 2:
                *(Uint16 *)p = pixel;
                break;
            case 3:
                if(SDL_BYTEORDER == SDL_BIG_ENDIAN) {
                    p[0] = (pixel >> 16) & 0xff;
                    p[1] = (pixel >> 8) & 0xff;
                    p[2] = pixel & 0xff;
                } else {
                    p[0] = pixel & 0xff;
                    p[1] = (pixel >> 8) & 0xff;
                    p[2] = (pixel >> 16) & 0xff;
                }
                break;
            case 4:
                *(Uint32 *)p = pixel;
                break;
        }
    }
/*----------------------------------------------------------------------------*/
int file_getFileName(int argc, char **argv, char *fileName)
{
    /* TODO: Implement */
    return 0;
}
/*----------------------------------------------------------------------------*/
int file_init(char *fileName, FILE **fileHandle, 
        int *noParticles, double *vector)
{
    /* TODO: Implement */
    *fileHandle = 0;
    return 0;
}
/*----------------------------------------------------------------------------*/
int file_close(FILE *fileHandle)
{
    if(fileHandle != 0)
    {
        close(fileHandle);
    }
    return 0;
}
/*----------------------------------------------------------------------------*/
int file_getNextParticleVector(FILE *fileHandle, int noParticles,
        double *vector)
{ 
    /* TODO: Implement */
    return 0;
}
/*------------------------------------------------------------------------------
 * Video Stuff
 *----------------------------------------------------------------------------*/
SDL_Surface * surface;
SDL_Surface * temp_surface;
/*----------------------------------------------------------------------------*/
int display_init()
{
    SDL_Init(SDL_INIT_VIDEO);
    SDL_WM_SetCaption("Gravisim Viewer", "Gravisim Viewer");
    surface = SDL_SetVideoMode(xResolution, yResolution, 0, 0);
    temp_surface = SDL_SetVideoMode(xResolution, yResolution, 0, 0);
    return 0;
}
/*----------------------------------------------------------------------------*/
int display_displayVector(int noParticles, double * vector)
{
    /* TODO: Implement */
    return 0;
}
/*----------------------------------------------------------------------------*/
int display_close()
{
    SDL_FreeSurface(surface);
    SDL_FreeSurface(temp_surface);
    SDL_Quit();
    return 0;
}
/*------------------------------------------------------------------------------
 * Control frames per second
 *----------------------------------------------------------------------------*/
struct timespec *timeToSleep;
struct timespec *timeRemaining;
int delay_init()
{
    timeToSleep = (struct timespec *)malloc(sizeof(struct timespec));
    timeRemaining = (struct timespec *)malloc(sizeof(struct timespec));
    timeToSleep->tv_sec = 0;
    timeToSleep->tv_nsec = 1000000000 / framesPerSecond;
    return 0;
}
/*----------------------------------------------------------------------------*/
int delay_close()
{
    free(timeToSleep);
    free(timeRemaining);
    timeToSleep = 0;
    timeRemaining = 0;
    return 0;
}
/*----------------------------------------------------------------------------*/
int delay()
{
    nanosleep(timeToSleep, timeRemaining);
    return 0;
}
/*------------------------------------------------------------------------------
 * MAIN
 *----------------------------------------------------------------------------*/
int main(int argc, char **argv)
{
    double  time;
    double *particleVector;
    char   *fileName;
    FILE   *fileHandle;
    int     noParticles;
    if( 0 != file_getFileName(argc, argv, fileName)) 
    {
        fprintf(stderr, "No file name given\n");
        exit(EXIT_FAILURE);
    }
    if(0 != file_init(fileName, &fileHandle, 
                &noParticles, particleVector)) 
    {
        fprintf(stderr, "Could not read in %s\n", fileName);
        exit(EXIT_FAILURE);
    }
    delay_init();
    display_init();
    do
    {
        display_displayVector(noParticles, particleVector);
        delay();
    } while(GOT_NEW_VECTOR == 
            file_getNextParticleVector(fileHandle, 
                noParticles, particleVector));
    file_close(fileHandle);
    /* free(particleVector); */
    display_close();
    exit(EXIT_SUCCESS);
}
/*
    if(SDL_Init(SDL_INIT_VIDEO) != 0) {
        cerr << "Konnte SDL nicht initialisieren!" << std::endl;
        return -1;
    };
    atexit(SDL_Quit);
    screen = SDL_SetVideoMode(640, 480, 16, 0);
    Uint32 white = SDL_MapRGB(screen->format, 0, 0, 0);
    SDL_LockSurface(screen);
    for(int i = 1; i < LEN; i++) {
        setPixel(screen, point2D[i][0], point2D[i][1], white);
    };
    SDL_UnlockSurface(screen);
    SDL_UpdateRect(screen, 0, 0, 0, 0);
    SDL_Delay(50000);
    cin >> c;
    free(screen);
    for( int i = 0; i < LEN; i++)
    {
        delete point2D[i];
    }
    delete[] point2D;
    delete trans;
    return 0;
}
*/
