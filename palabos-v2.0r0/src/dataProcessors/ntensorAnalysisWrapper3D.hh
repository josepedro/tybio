/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file
 * Helper functions for domain initialization -- header file.
 */
#ifndef NTENSOR_ANALYSIS_WRAPPER_3D_HH
#define NTENSOR_ANALYSIS_WRAPPER_3D_HH

#include "dataProcessors/ntensorAnalysisWrapper3D.h"
#include "dataProcessors/ntensorAnalysisFunctional3D.h"
#include "atomicBlock/reductiveDataProcessorWrapper3D.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"
#include "multiBlock/multiBlockGenerator3D.h"

namespace plb {

template<typename T1, typename T2>
void copy( MultiNTensorField3D<T1>& field,
           MultiNTensorField3D<T2>& convertedField, Box3D domain)
{
    applyProcessingFunctional (
            new CopyConvertNTensorFunctional3D<T1,T2>, domain, field, convertedField );
}

template<typename T1, typename T2>
MultiNTensorField3D<T2>* copyConvert( MultiNTensorField3D<T1>& field,
                                      Box3D domain)
{
    MultiNTensorField3D<T2>* convertedField
        = generateMultiNTensorField<T2>(field, domain, field.getNdim());
    plb::copy(field, *convertedField, domain);
    return convertedField;
}

template<typename T>
void nTensorToScalar(MultiNTensorField3D<T>& nTensor, MultiScalarField3D<T>& scalar) {
    PLB_ASSERT( nTensor.getNdim()==1 );
    std::vector<MultiBlock3D*> args;
    args.push_back(&nTensor);
    args.push_back(&scalar);
    applyProcessingFunctional (
            new ConvertNTensorToScalarFunctional3D<T>, nTensor.getBoundingBox(), args );
}

template<typename T>
void scalarToNTensor(MultiScalarField3D<T>& scalar, MultiNTensorField3D<T>& nTensor) {
    std::vector<MultiBlock3D*> args;
    args.push_back(&scalar);
    args.push_back(&nTensor);
    applyProcessingFunctional (
            new ConvertScalarToNTensorFunctional3D<T>, scalar.getBoundingBox(), args );
}

template<typename T>
std::auto_ptr<MultiNTensorField3D<T> > scalarToNTensor(MultiScalarField3D<T>& scalar) {
    plint nDim = 1;
    std::auto_ptr<MultiNTensorField3D<T> > nTensor =
        defaultGenerateMultiNTensorField3D<T>(scalar.getMultiBlockManagement(), nDim);
    scalarToNTensor(scalar, *nTensor);
    return nTensor;
}

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > nTensorToScalar(MultiNTensorField3D<T>& nTensor) {
    std::auto_ptr<MultiScalarField3D<T> > scalar =
        generateMultiScalarField<T>(nTensor, nTensor.getBoundingBox());
    nTensorToScalar(nTensor, *scalar);
    return scalar;
}

}  // namespace plb

#endif  // NTENSOR_ANALYSIS_WRAPPER_3D_HH
