//
// Billion Star 3D Rendering Engine
// Kevin M. Loch
//
// 3D rendering engine for the ESA Gaia EDR3 star dataset

/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2021, Kevin Loch
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "bsrender.h" // needs to be first to get GNU_SOURCE define for strcasestr
#include <stdio.h>
#include <stdlib.h>
#include <sys/wait.h>

int checkExceptions(bsr_state_t *bsr_state) {
  int i;

  if (bsr_state->perthread->my_thread_id == 0) {
    // main thread

    // TBD: see if httpd parent's tcp socket is in CLOSE_WAIT
    // This could be caused by closing browser/tab or hitting stop button

    // see if parent httpd process has died
    if (getppid() != bsr_state->httpd_pid) {
      // parent httpd process has died, exit
      exit(1);
    }
    
    // clean up zombie processes
    waitpid(-1, NULL, WNOHANG);

    // check for child processes that have died
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      if (getpgid(bsr_state->status_array[i].pid) != bsr_state->master_pgid) {
        // if any child process has died we can't finish, exit
        exit(1);
      }
    }
  } else {
    // worker thread

    // check if main thread has died
    if (getppid() != bsr_state->master_pid) {
      // main thread has died, exit
      exit(1);
    }
  }

  return(0);
}

int waitForWorkerThreads(bsr_state_t *bsr_state, int min_status) {
  int i;
  int loop_count;
  volatile int all_workers_done;

  loop_count=0;
  all_workers_done=0;
  while (all_workers_done == 0) {

    // periodically check for exceptions
    loop_count++;
    if ((loop_count % 10000) == 0) {
      checkExceptions(bsr_state);
      loop_count=1;
    }

    // see if all worker threads have completed task
    all_workers_done=1;
    for (i=1; i <= bsr_state->num_worker_threads; i++) {
      if (bsr_state->status_array[i].status < min_status) {
        all_workers_done=0;
      }
    }
  }

  return(0);
}

int waitForMainThread(bsr_state_t *bsr_state, int min_status) {
  volatile int cont;
  int loop_count;

  loop_count=0;
  cont=0;
  while (cont == 0) {

    // periodically check for exceptions
    loop_count++;
    if ((loop_count % 10000) == 0) {
      checkExceptions(bsr_state);
      loop_count=1;
    }

    // see if main thread says go to next task
    if (bsr_state->status_array[bsr_state->perthread->my_thread_id].status >= min_status) {
      cont=1;
    }
  }

  return(0);
}

int limitIntensity(double *pixel_r, double *pixel_g, double *pixel_b) {

  //
  // limit pixel to range 0.0-1.0 without regard to color
  //
  if (*pixel_r < 0.0) {
    *pixel_r=0.0;
  }
  if (*pixel_g < 0.0) {
    *pixel_g=0.0;
  }
  if (*pixel_b < 0.0) {
    *pixel_b=0.0;
  }

  if (*pixel_r > 1.0) {
    *pixel_r=1.0;
  }
  if (*pixel_g > 1.0) {
    *pixel_g=1.0;
  }
  if (*pixel_b > 1.0) {
    *pixel_b=1.0;
  }

  return(0);
}

int limitIntensityPreserveColor(double *pixel_r, double *pixel_g, double *pixel_b) {
  double pixel_max;

  //
  // limit pixel to range 0.0-1.0 while maintaining color (max channel=1.0)
  //
  if (*pixel_r < 0.0) {
    *pixel_r=0;
  }
  if (*pixel_g < 0.0) {
    *pixel_g=0;
  }
  if (*pixel_b < 0.0) {
    *pixel_b=0;
  }

  if ((*pixel_r > 1.0) || (*pixel_g > 1.0) || (*pixel_b > 1.0)) {
    pixel_max=0;
    if (*pixel_r > pixel_max) {
      pixel_max=*pixel_r;
    }
    if (*pixel_g > pixel_max) {
      pixel_max=*pixel_g;
    }
    if (*pixel_b > pixel_max) {
      pixel_max=*pixel_b;
    }
    *pixel_r=*pixel_r / pixel_max;
    *pixel_g=*pixel_g / pixel_max;
    *pixel_b=*pixel_b / pixel_max;
  }

  return(0);
}
