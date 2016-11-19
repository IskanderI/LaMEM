/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **    Copyright (c) 2011-2015, JGU Mainz, Anton Popov, Boris Kaus
 **    All rights reserved.
 **
 **    This software was developed at:
 **
 **         Institute of Geosciences
 **         Johannes-Gutenberg University, Mainz
 **         Johann-Joachim-Becherweg 21
 **         55128 Mainz, Germany
 **
 **    project:    LaMEM
 **    filename:   marker.h
 **
 **    LaMEM is free software: you can redistribute it and/or modify
 **    it under the terms of the GNU General Public License as published
 **    by the Free Software Foundation, version 3 of the License.
 **
 **    LaMEM is distributed in the hope that it will be useful,
 **    but WITHOUT ANY WARRANTY; without even the implied warranty of
 **    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 **    See the GNU General Public License for more details.
 **
 **    You should have received a copy of the GNU General Public License
 **    along with LaMEM. If not, see <http://www.gnu.org/licenses/>.
 **
 **
 **    Contact:
 **        Boris Kaus       [kaus@uni-mainz.de]
 **        Anton Popov      [popov@uni-mainz.de]
 **
 **
 **    Main development team:
 **         Anton Popov      [popov@uni-mainz.de]
 **         Boris Kaus       [kaus@uni-mainz.de]
 **         Tobias Baumann
 **         Adina Pusok
 **         Arthur Bauville
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
//---------------------------------------------------------------------------
//........................   MARKER ROUTINES   ..............................
//---------------------------------------------------------------------------
#ifndef __marker_h__
#define __marker_h__
//---------------------------------------------------------------------------

// input polygon data
typedef struct
{
	PetscInt     dir; // normal vector of polygon plane
	PetscInt   ax[2]; // axis that span the polygon plane
	PetscInt   phase; // phase that the polygon defines
	PetscInt    type; // type can be of additive or assigning nature
	PetscInt     num; // number of polygon slices defining the volume
	PetscInt     len; // number of nodes of polygon
	PetscInt    idxs; // index of first polygon slice
	PetscInt    gidx; // global plane index (consistent with planes of markers)
	PetscInt    lidx; // local plane index (consistent with planes of markers)
	PetscInt       n; // number of polygon nodes
	PetscInt   nmark; // number of markers in current volume
	PetscScalar   *X; // coordinates of the polygon x1,y1,x2,y2,...xn,yn

} Polygon2D;

//---------------------------------------------------------------------------

// markers initialization
PetscErrorCode ADVMarkInit(AdvCtx *actx);

// generate coordinates of uniformly distributed markers
PetscErrorCode ADVMarkInitCoord(AdvCtx *actx);

// save all local markers to disk (parallel output)
PetscErrorCode ADVMarkSave(AdvCtx *actx);

// check phase IDs of all the markers
PetscErrorCode ADVMarkCheckMarkers(AdvCtx *actx);

//---------------------------------------------------------------------------

// Specific initialization routines

PetscErrorCode ADVMarkInitFileParallel (AdvCtx *actx);
PetscErrorCode ADVMarkInitFileRedundant(AdvCtx *actx);
PetscErrorCode ADVMarkInitFilePolygons (AdvCtx *actx);
PetscErrorCode ADVMarkInitDiapir       (AdvCtx *actx);
PetscErrorCode ADVMarkInitBlock        (AdvCtx *actx);
PetscErrorCode ADVMarkInitSubduction   (AdvCtx *actx);
PetscErrorCode ADVMarkInitFolding      (AdvCtx *actx);
PetscErrorCode ADVMarkInitDetachment   (AdvCtx *actx);
PetscErrorCode ADVMarkInitSlab         (AdvCtx *actx);
PetscErrorCode ADVMarkInitSpheres      (AdvCtx *actx);
PetscErrorCode ADVMarkInitBands        (AdvCtx *actx);
PetscErrorCode ADVMarkInitDomes        (AdvCtx *actx);
PetscErrorCode ADVMarkInitRotation     (AdvCtx *actx);

//---------------------------------------------------------------------------

// service functions

PetscErrorCode ADVMarkSetTempFromFile  (AdvCtx *actx);

void ADVMarkSecIdx(AdvCtx *actx, PetscInt dir, PetscInt Nslice, PetscInt *idx);

//---------------------------------------------------------------------------

// definitions

#ifndef max
    #define max(a,b) (a >= b ? a : b)
    #define min(a,b) (a <= b ? a : b)
#endif

//---------------------------------------------------------------------------
#endif
