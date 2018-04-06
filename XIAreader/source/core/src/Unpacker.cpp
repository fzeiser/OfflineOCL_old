/*******************************************************************************
 * Copyright (C) 2016 Vetle W. Ingeberg                                        *
 * Author: Vetle Wegner Ingeberg, v.w.ingeberg@fys.uio.no                      *
 *                                                                             *
 * --------------------------------------------------------------------------- *
 * This program is free software; you can redistribute it and/or modify it     *
 * under the terms of the GNU General Public License as published by the       *
 * Free Software Foundation; either version 3 of the license, or (at your      *
 * option) any later version.                                                  *
 *                                                                             *
 * This program is distributed in the hope that it will be useful, but         *
 * WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General   *
 * Public License for more details.                                            *
 *                                                                             *
 * You should have recived a copy of the GNU General Public License along with *
 * the program. If not, see <http://www.gnu.org/licenses/>.                    *
 *                                                                             *
 *******************************************************************************/

/*!
 * \file Unpacker.cpp
 * \brief Implementation of Unpacker.
 * \author Vetle W. Ingeberg
 * \version 0.8.0
 * \date 2015-2016
 * \copyright GNU Public License v. 3
 */

#include "Unpacker.h"

#include "TDRWordBuffer.h"
#include "Event.h"
#include "experimentsetup.h"
#include "DefineFile.h"

#include <iostream>
#include <string>
#include <cstdlib>

#define GAP_SIZE 100

Unpacker::Unpacker()
	: buffer( 0 )
	, buffer_idx( 0 )
	, eventlength_sum( 0 )
	, event_count ( 0 )
{
}

void Unpacker::SetBuffer(const WordBuffer* buffr)
{
    buffer = buffr; // Setting internal variable with the new buffer.


    buffer_idx = 0; // Reset buffer position.
    curr_Buf = 0; // Reseting the position in the buffer.

    // Resting the length of the buffer.
	event_count = eventlength_sum = 0;
}

Unpacker::Status Unpacker::Next(Event &event)
{
    if ( curr_Buf >= buffer->GetSize() )
        return END; // End of buffer reached.

	int n_data = 0;
    if ( !UnpackOneEvent(event, n_data) ){
        buffer_idx = n_data + 1;
        return END; // If no event found, end of buffer.
	}

    eventlength_sum += event.length;
	event_count += 1;
    buffer_idx = n_data + 1;

	return OKAY;
}
#if SINGLES
bool Unpacker::UnpackOneEvent(Event& event, int& n_data)
{

    event.Reset();

        if ( curr_Buf >= buffer->GetSize() )
            return false;

        int64_t timediff;
        int start = curr_Buf;
        int stop = curr_Buf+1;
        subevent_t sevt1 = subevent_t((*buffer)[curr_Buf]), sevt2;
        for (size_t i = curr_Buf + 1 ; i < buffer->GetSize() ; ++i){
            sevt2 = sevt1;
            sevt1 = subevent_t((*buffer)[i]);
            timediff = sevt1.timestamp - sevt2.timestamp;
            if (timediff > GAP_SIZE){
                stop = i;
                curr_Buf = stop;
                return PackEvent(event, start, stop);
            }
        }

        stop = buffer->GetSize();
        curr_Buf = stop;
        return PackEvent(event, start, stop);
}
#else
bool Unpacker::UnpackOneEvent(Event& event, int& n_data)
{
    event.Reset();

    if ( curr_Buf >= buffer->GetSize() )
        return false;

    int64_t timediff;
    word_t cWord;
    for (int i =  curr_Buf ; i < buffer->GetSize() ; ++i){
        cWord = (*buffer)[i];
        if (GetDetector(cWord.address).type == ppac){
            int start = i;
            int stop = i+1;

            // We will move back a given number of words and check timestamps
            for (int j = i ; j > 0 ; --j){
                timediff = (*buffer)[j-1].timestamp - cWord.timestamp;
                if (abs(timediff) > 2500){
                    start = j;
                    break;
                }
            }

            for (int j = i ; j < buffer->GetSize() - 1 ; ++j){
                timediff = (*buffer)[j+1].timestamp - cWord.timestamp;
                if (abs(timediff) > 2500){
                    stop = j+1;
                    break;
                }
            }

            curr_Buf = stop;
            n_data = stop-start;
            return PackEvent(event, start, stop);
        }
    }

    curr_Buf = buffer->GetSize();
    return false;
}
#endif // SINGLES

bool Unpacker::PackEvent(Event& event, int start, int stop)
{
    event.length = stop - start;
    DetectorInfo_t dinfo;
    for (int i = start ; i < stop ; ++i){
        dinfo = GetDetector((*buffer)[i].address);

        switch (dinfo.type) {
        case labr: {
            if ( event.n_labr[dinfo.detectorNum] < MAX_WORDS_PER_DET &&
                 dinfo.detectorNum < NUM_LABR_DETECTORS){
                event.w_labr[dinfo.detectorNum][event.n_labr[dinfo.detectorNum]++] = subevent_t((*buffer)[i]);
                ++event.tot_labr;
            } else {
                std::cerr << __PRETTY_FUNCTION__ << ": Could not populate LaBr word, run debugger with appropriate break point for more details" << std::endl;
            }
            break;
        }
        case deDet: {
            if ( event.n_dEdet[dinfo.detectorNum] < MAX_WORDS_PER_DET &&
                 dinfo.detectorNum < NUM_SI_DE_DET){
                event.w_dEdet[dinfo.detectorNum][event.n_dEdet[dinfo.detectorNum]++] = subevent_t((*buffer)[i]);
                ++event.tot_dEdet;
            } else {
                std::cerr << __PRETTY_FUNCTION__ << ": Could not populate dEdet word, run debugger with appropriate break point for more details" << std::endl;
            }
            break;
        }
        case eDet: {
            if ( event.n_Edet[dinfo.detectorNum] < MAX_WORDS_PER_DET &&
                 dinfo.detectorNum < NUM_SI_E_DET){
                event.w_Edet[dinfo.detectorNum][event.n_Edet[dinfo.detectorNum]++] = subevent_t((*buffer)[i]);
                ++event.tot_Edet;
            } else {
                std::cerr << __PRETTY_FUNCTION__ << ": Could not populate Edet word, run debugger with appropriate break point for more details" << std::endl;
            }
            break;
        }
        case eGuard: {
            if ( event.n_Eguard[dinfo.detectorNum] < MAX_WORDS_PER_DET &&
                 dinfo.detectorNum < NUM_SI_E_GUARD){
                event.w_Eguard[dinfo.detectorNum][event.n_Eguard[dinfo.detectorNum]++] = subevent_t((*buffer)[i]);
                ++event.tot_Eguard;
            } else {
                std::cerr << __PRETTY_FUNCTION__ << ": Could not populate eGuard word, run debugger with appropriate break point for more details" << std::endl;
            }
            break;
        }
        case ppac: {
            if ( event.n_ppac[dinfo.detectorNum] < MAX_WORDS_PER_DET &&
                 dinfo.detectorNum < NUM_PPAC){
                event.w_ppac[dinfo.detectorNum][event.n_ppac[dinfo.detectorNum]++] = subevent_t((*buffer)[i]);
                ++event.tot_ppac;
            } else {
                std::cerr << __PRETTY_FUNCTION__ << ": Could not populate PPAC word, run debugger with appropriate break point for more details" << std::endl;
            }
            break;
        }
        case rfchan: {
            if ( event.n_RFpulse < MAX_WORDS_PER_DET )
                event.w_RFpulse[event.n_RFpulse++] = subevent_t((*buffer)[i]);
            else
                std::cerr << __PRETTY_FUNCTION__ << ": Could not populate RF word, run debugger with appropriate break point for more details" << std::endl;
            break;
        }
        default:
            break;
        }
    }
    return true;
}
