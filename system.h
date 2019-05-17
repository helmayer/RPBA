/*

RPBA - Robust Parallel Bundle Adjustment

File system.h



Copyright 2019 Mario Michelini, Bundeswehr University Munich, Germany, Mario.Michelini@unibw.de

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/


#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#ifdef _WIN32
	#include <windows.h>
#else
	#include <sys/time.h>
	#include <unistd.h>
#endif
#ifdef _OPENMP
	#include <omp.h>
#endif
#include <thread>


/**
 * System interface
 * @author Mario Michelini
 */
class System
{
public:
	typedef std::size_t Time;

public:
	/**
	 * Returns the maximum number of threads
	 * @return Number of threads
	 */
	static int threadCount()
	{
		#ifdef _OPENMP
			int mt = omp_get_num_procs();
			if( mt > 0 )
				return mt;
		#endif

		mt = (int) std::thread::hardware_concurrency();
		if( mt > 0 )
			return mt;

		#ifdef _WIN32
			const char* evar = getenv("NUMBER_OF_PROCESSORS");
			if( evar )
			{
				mt = std::stoi(evar);
				if( mt > 0 )
					return mt;
			}
		#else
			mt = (int) sysconf(_SC_NPROCESSORS_ONLN);
			if( mt > 0 )
				return mt;
		#endif

		return 1;
	}

	/**
	 * Retrieves a time-stamp in milliseconds.
	 * The difference between two time-stamps specifies the elapsed time in milliseconds.
	 * @return Time-stamp in milliseconds
	 */
	static Time getTickCount()
	{
		#ifdef _WIN32
			return (Time) GetTickCount();
		#else
			timeval ts;
			gettimeofday(&ts,0);
			return (Time) ts.tv_sec * 1000 + (ts.tv_usec / 1000);
		#endif
	}

	/** Timer */
	class Timer
	{
		System::Time m_accum;	/**< Accumulated time */
		System::Time m_start;	/**< Start time-stamp */

	public:
		Timer() { reset(); }

		/** Starts the timer */
		inline void start() { m_start = System::getTickCount(); }

		/** Stops the timer */
		inline void stop()
		{
			m_accum += System::getTickCount() - m_start;
			m_start = 0;
		}

		/** Resets the timer */
		void reset()
		{
			m_accum = 0;
			start();
		}
	};
};


#endif // _SYSTEM_H_
