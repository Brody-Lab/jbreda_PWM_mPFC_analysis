# jbreda_PWM_ephys_analysis

repository for analysis of PWM ephys data (as of 2020_1_6 this is for wireless mPFC tetrode recordings)

# TODO
* change naming to match tyler;s
* add plots to find_wireless_sess_J
* rerun all sessions with scopeon alignment

## Assumptions
1. This repo is built as a continuation of the output from the [jbreda_kilsort](https://github.com/Brody-Lab/jbreda_kilosort) repository that was designed for analyzing wireless tetrode data in kilosort. This is the best place to look if you want to understand how my files are structured. Additionally, I have taken & adjusted function from Tyler and [Emily's](https://github.com/Brody-Lab/emilyjanedennis_PWManalysis/blob/master/find_wireless_sess.m) analyses repositories.

2. While you could implement some parts of this code pre or during sorting, I am implementing it after all my sorting has been completed for one rat. So, all the sessions I am interested in (from raw, to preprocessed, to sorted) a convienently represented together in their respected directories. Thereby, most of these functions are intended to run *all* sessions for *one* rat as opposed to *some* sessions across *many* rats. Note that adjusting this should be quite easy.


## Session alignment

1. `find_wireless_sess_J` takes a session and finds the corresponding behavioral session in bdata.

**specific details**:

For dealing with multiple sessions on the same day, this function compares TTL ITI from the headstage (trodes) and computer (finite state machine) and finds the session with the smallest residual after a linear regression. Additionally, using the output from this regression of headstage ttl vs computer ttls, one can now move from spike times to behavior times.

**output**

    Creates a structure with:
