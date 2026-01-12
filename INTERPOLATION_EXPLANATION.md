# Interpolation Weights Explanation and Marker-to-Face Interpolation

## How the Weights are Calculated (Lines 2078-2085)

### Weight Definitions

The weights are calculated using two macros defined in `src/fdstag.h`:

1. **WEIGHT_POINT_CELL(i, x, ds)**: 
   ```cpp
   (1.0 - PetscAbsScalar(x - ds.ccoor[i])/(ds.ncoor[i+1] - ds.ncoor[i]))
   ```
   - Calculates the interpolation weight for a point `x` within cell `i`'s control volume
   - Returns 1.0 at the cell center (`ds.ccoor[i]`), decreasing linearly to 0 at cell boundaries
   - Denominator: cell size = distance between bounding nodes (`ds.ncoor[i+1] - ds.ncoor[i]`)
   - This is a **linear distance-based weight** within the cell

2. **WEIGHT_POINT_NODE(i, x, ds)**:
   ```cpp
   (1.0 - PetscAbsScalar(x - ds.ncoor[i])/(ds.ccoor[i] - ds.ccoor[i-1]))
   ```
   - Calculates the interpolation weight for a point `x` within node `i`'s control volume
   - Returns 1.0 at the node (`ds.ncoor[i]`), decreasing linearly to 0 at control volume boundaries
   - Denominator: node control volume size = distance between neighboring cell centers
   - This is a **linear distance-based weight** within the node's control volume

### What These Weights Are Doing

In the code at lines 2078-2085:

1. **Cell weights (wxc, wyc, wzc)**: Measure how close the marker is to the cell center in each direction
   - Used for interpolating to cell-centered quantities
   - Range: [0, 1], where 1 = at cell center, 0 = at cell boundary

2. **Node weights (wxn, wyn, wzn)**: Measure how close the marker is to the edge node in each direction
   - Used for interpolating to edge nodes (DA_XY, DA_XZ, DA_YZ)
   - Range: [0, 1], where 1 = at node, 0 = at control volume boundary

3. **Combined weights** (lines 2091-2093):
   - `wxn*wyn*wzc` for XY edges: X and Y node weights, Z cell weight
   - `wxn*wyc*wzn` for XZ edges: X node, Y cell, Z node weights  
   - `wxc*wyn*wzn` for YZ edges: X cell, Y and Z node weights
   - These combinations distribute marker data to the appropriate edge locations

## Grid Structure

The FDSTAG grid uses a **staggered grid** layout:
- **DA_CEN**: Cell centers (pressure, temperature)
- **DA_X, DA_Y, DA_Z**: Face centers (velocities normal to faces)
- **DA_XY, DA_XZ, DA_YZ**: Edge nodes (shear stresses)
- **DA_COR**: Corner nodes

Face locations:
- **DA_X faces**: Located at X-nodes, Y-cell centers, Z-cell centers
- **DA_Y faces**: Located at X-cell centers, Y-nodes, Z-cell centers  
- **DA_Z faces**: Located at X-cell centers, Y-cell centers, Z-nodes

## Interpolating from Markers to Faces

Since markers are a cloud of points and you get the closest cell with `get_cell_ijk`, here are **three options** for interpolating marker data to faces:

---

## Option 1: Direct Weight-Based Interpolation (Similar to Edge Interpolation)

This approach uses the same weight calculation philosophy but adapts it for faces.

### For DA_X faces (X-normal faces):
```cpp
// After getting cell I, J, K from get_cell_ijk
// Get marker coordinates
xp = P->X[0];
yp = P->X[1];
zp = P->X[2];

// Get cell center
xc = fs->dsx.ccoor[I];
yc = fs->dsy.ccoor[J];
zc = fs->dsz.ccoor[K];

// Find which X-face the marker is closest to
if(xp > xc) { II = I+1; } else { II = I; }

// Calculate weights
wxn = WEIGHT_POINT_NODE(II, xp, fs->dsx);  // X-direction: node weight
wyc = WEIGHT_POINT_CELL(J, yp, fs->dsy);   // Y-direction: cell weight
wzc = WEIGHT_POINT_CELL(K, zp, fs->dsz);  // Z-direction: cell weight

// Interpolate to X-face
lvx[sz+K][sy+J][sx+II] += wxn*wyc*wzc*P->value;
```

### For DA_Y faces (Y-normal faces):
```cpp
// Find which Y-face
if(yp > yc) { JJ = J+1; } else { JJ = J; }

wxc = WEIGHT_POINT_CELL(I, xp, fs->dsx);   // X-direction: cell weight
wyn = WEIGHT_POINT_NODE(JJ, yp, fs->dsy);  // Y-direction: node weight
wzc = WEIGHT_POINT_CELL(K, zp, fs->dsz);   // Z-direction: cell weight

lvy[sz+K][sy+JJ][sx+I] += wxc*wyn*wzc*P->value;
```

### For DA_Z faces (Z-normal faces):
```cpp
// Find which Z-face
if(zp > zc) { KK = K+1; } else { KK = K; }

wxc = WEIGHT_POINT_CELL(I, xp, fs->dsx);   // X-direction: cell weight
wyc = WEIGHT_POINT_CELL(J, yp, fs->dsy);   // Y-direction: cell weight
wzn = WEIGHT_POINT_NODE(KK, zp, fs->dsz);  // Z-direction: node weight

lvz[sz+KK][sy+J][sx+I] += wxc*wyc*wzn*P->value;
```

**Pros**: Simple, consistent with existing edge interpolation
**Cons**: Only uses nearest face, may not be as accurate for markers far from faces

---

## Option 2: Trilinear Interpolation (Like InterpLin3D)

This uses the same approach as `InterpLin3D` but in reverse (marker → grid).

### For DA_X faces:
Each X-face is surrounded by 4 X-faces (2 in Y, 2 in Z directions). You need to:
1. Find the 4 surrounding X-faces
2. Calculate trilinear weights
3. Distribute marker value to all 4 faces

```cpp
// Get marker coordinates
xp = P->X[0];
yp = P->X[1];
zp = P->X[2];

// Get cell center
xc = fs->dsx.ccoor[I];
yc = fs->dsy.ccoor[J];
zc = fs->dsz.ccoor[K];

// Find base face indices
if(xp > xc) { II = I+1; } else { II = I; }

// Get coordinates of surrounding faces
PetscScalar *ncx = fs->dsx.ncoor;  // X-face locations
PetscScalar *ccy = fs->dsy.ccoor;  // Y-cell centers
PetscScalar *ccz = fs->dsz.ccoor;  // Z-cell centers

// Calculate relative coordinates within the interpolation cube
PetscScalar xe = (xp - ncx[II])/(ncx[II+1] - ncx[II]); 
PetscScalar xb = 1.0 - xe;
PetscScalar ye = (yp - ccy[J])/(ccy[J+1] - ccy[J]);
PetscScalar yb = 1.0 - ye;
PetscScalar ze = (zp - ccz[K])/(ccz[K+1] - ccz[K]);
PetscScalar zb = 1.0 - ze;

// Distribute to 4 surrounding X-faces (2x2 in Y-Z plane)
lvx[sz+K  ][sy+J  ][sx+II] += xb*yb*zb*P->value;
lvx[sz+K  ][sy+J  ][sx+II+1] += xe*yb*zb*P->value;  // if II+1 valid
lvx[sz+K  ][sy+J+1][sx+II] += xb*ye*zb*P->value;
lvx[sz+K  ][sy+J+1][sx+II+1] += xe*ye*zb*P->value;  // if II+1 valid
lvx[sz+K+1][sy+J  ][sx+II] += xb*yb*ze*P->value;
lvx[sz+K+1][sy+J  ][sx+II+1] += xe*yb*ze*P->value;  // if II+1 valid
lvx[sz+K+1][sy+J+1][sx+II] += xb*ye*ze*P->value;
lvx[sz+K+1][sy+J+1][sx+II+1] += xe*ye*ze*P->value;  // if II+1 valid
```

**Note**: For X-faces, X varies along nodes, Y and Z vary along cell centers. Adjust indices accordingly.

**Pros**: More accurate, smooth interpolation
**Cons**: More complex, requires careful index management

---

## Option 3: Inverse Distance Weighting from Multiple Markers

Instead of interpolating one marker at a time, collect all markers near each face and use inverse distance weighting.

```cpp
// For each face, find nearby markers and weight by inverse distance
for each face (i, j, k) {
    PetscScalar sum_weight = 0.0;
    PetscScalar sum_value = 0.0;
    
    for each marker P in nearby cells {
        PetscScalar dist = distance(face_location, P->X);
        PetscScalar weight = 1.0 / (dist + epsilon);  // epsilon prevents division by zero
        sum_weight += weight;
        sum_value += weight * P->value;
    }
    
    face_value[i][j][k] = sum_value / sum_weight;
}
```

**Pros**: Handles sparse marker distributions well
**Cons**: More computationally expensive, requires neighbor search

---

## Recommended Approach

For consistency with the existing codebase, I recommend **Option 1** (Direct Weight-Based) because:
1. It's consistent with the existing edge interpolation pattern
2. It's simple and efficient
3. It uses the same weight calculation macros
4. It follows the same logic as lines 2078-2093

However, if you need higher accuracy, **Option 2** (Trilinear) is better but requires more careful implementation.

---

## Questions to Clarify

1. **What field are you interpolating?** (velocity, stress, temperature, etc.)
2. **Do you need to accumulate from multiple markers?** (like the `+=` in line 2091)
3. **Do you need normalization?** (dividing by total weight after accumulation)
4. **Which faces specifically?** (DA_X, DA_Y, DA_Z, or all three?)
5. **Is this for a specific function or general-purpose interpolation?**

---

## Example Implementation Template

Here's a template function following Option 1:

```cpp
PetscErrorCode ADVInterpMarkToFace(AdvCtx *actx, PetscInt face_type)
{
    FDSTAG *fs = actx->fs;
    JacRes  *jr = actx->jr;
    Marker  *P;
    
    PetscInt     jj, ID, I, J, K, II, JJ, KK;
    PetscInt     sx, sy, sz, nx, ny;
    PetscScalar  xp, yp, zp, xc, yc, zc;
    PetscScalar  wxc, wyc, wzc, wxn, wyn, wzn;
    PetscScalar ***lvx, ***lvy, ***lvz;
    
    // Get starting indices
    sx = fs->dsx.pstart; nx = fs->dsx.ncels;
    sy = fs->dsy.pstart; ny = fs->dsy.ncels;
    sz = fs->dsz.pstart;
    
    // Access face arrays
    ierr = DMDAVecGetArray(fs->DA_X, jr->lvx, &lvx); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(fs->DA_Y, jr->lvy, &lvy); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(fs->DA_Z, jr->lvz, &lvz); CHKERRQ(ierr);
    
    // Zero the arrays if accumulating
    // ierr = VecZeroEntries(jr->lvx); CHKERRQ(ierr);
    // ... etc
    
    // Loop over markers
    for(jj = 0; jj < actx->nummark; jj++)
    {
        P = &actx->markers[jj];
        
        // Get host cell
        ID = actx->cellnum[jj];
        GET_CELL_IJK(ID, I, J, K, nx, ny)
        
        // Get marker coordinates
        xp = P->X[0];
        yp = P->X[1];
        zp = P->X[2];
        
        // Get cell center
        xc = fs->dsx.ccoor[I];
        yc = fs->dsy.ccoor[J];
        zc = fs->dsz.ccoor[K];
        
        // Interpolate to X-faces
        if(face_type == 0 || face_type == -1) {  // -1 = all faces
            if(xp > xc) { II = I+1; } else { II = I; }
            wxn = WEIGHT_POINT_NODE(II, xp, fs->dsx);
            wyc = WEIGHT_POINT_CELL(J, yp, fs->dsy);
            wzc = WEIGHT_POINT_CELL(K, zp, fs->dsz);
            lvx[sz+K][sy+J][sx+II] += wxn*wyc*wzc*P->value;  // Replace with your field
        }
        
        // Interpolate to Y-faces
        if(face_type == 1 || face_type == -1) {
            if(yp > yc) { JJ = J+1; } else { JJ = J; }
            wxc = WEIGHT_POINT_CELL(I, xp, fs->dsx);
            wyn = WEIGHT_POINT_NODE(JJ, yp, fs->dsy);
            wzc = WEIGHT_POINT_CELL(K, zp, fs->dsz);
            lvy[sz+K][sy+JJ][sx+I] += wxc*wyn*wzc*P->value;
        }
        
        // Interpolate to Z-faces
        if(face_type == 2 || face_type == -1) {
            if(zp > zc) { KK = K+1; } else { KK = K; }
            wxc = WEIGHT_POINT_CELL(I, xp, fs->dsx);
            wyc = WEIGHT_POINT_CELL(J, yp, fs->dsy);
            wzn = WEIGHT_POINT_NODE(KK, zp, fs->dsz);
            lvz[sz+KK][sy+J][sx+I] += wxc*wyc*wzn*P->value;
        }
    }
    
    // Restore access
    ierr = DMDAVecRestoreArray(fs->DA_X, jr->lvx, &lvx); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(fs->DA_Y, jr->lvy, &lvy); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(fs->DA_Z, jr->lvz, &lvz); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}
```
