# rtSOUND plugin for NeuroSimo real-time software

Real-time EEG noise removal plugin for NeuroSimo real-time software. See the [NeuroSimo repository](https://github.com/NeuroSimo/neurosimo) for information and installation instructions.

## Getting started

1. [Create a project](https://github.com/NeuroSimo/neurosimo/blob/main/md/getting-started.md) in NeuroSimo.
2. Place `rtsound.py` and `leadfield.csv` in the `preprocessor` in the preprocessor folder.
3. Change the `leadfield.csv` file to a lead-field matrix compatible with your real-time-streaming EEG data.
4. Enable the preprocessor and select `rtsound` in the preprocessor dropdown menu.
5. If you wish to change rtSOUND parameters, they can be changed in the `rtsound.py` script section "User-changeable parameters".

## License
This project is licensed under the GPL v3 License - see the [LICENSE](https://github.com/matildamakkonen/rt-sound-neurosimo-plugin/blob/main/LICENSE) file for details.

## References
Mutanen, T. P., Metsomaa, J., Makkonen, M., Varone, G., Marzetti, L., & Ilmoniemi, R. J. (2022). Source-based artifact-rejection techniques for TMS–EEG. Journal of Neuroscience Methods, 382, 109693.

Makkonen, M., Mutanen, T., Metsomaa, J., Zrenner, C., Souza, V., & Ilmoniemi, R. (2021). Real-time artifact detection and removal for closed-loop EEG-TMS. International Journal of Bioelectromagnetism, 23(2), 1-4.

Mutanen, T. P., Metsomaa, J., Liljander, S., & Ilmoniemi, R. J. (2018). Automatic and robust noise suppression in EEG and MEG: The SOUND algorithm. Neuroimage, 166, 135-151.
