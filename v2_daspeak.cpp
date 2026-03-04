/*
 * DaSpeak Formant Synthesizer - C++ Implementation
 * Copyright (c) 2026 hotdogdevourer
 *
 * Licensed under the MIT License.
 * See the LICENSE file in the project root for full license information.
 */
 
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <cmath>
#include <random>
#include <fstream>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <iomanip>
constexpr double M_PI = 3.14159265358979323846;

namespace DaSpeakSynth {
    
const std::set<std::string> VOWELS = {
    "AH", "AE", "AA", "AO", "EH", "EY", "IH", "IY", "OW", "UH", "UW", "ER"
};

const std::set<std::string> STOPS = {
    "P", "T", "K", "B", "D", "G", "CH"
};

const std::set<std::string> VALID_PHONEMES = {
    "SIL", "AH", "AE", "AA", "AO", "EH", "EY", "IH", "IY", "OW", "UH", "UW", "ER",
    "B", "D", "G", "P", "T", "K", "M", "N", "NG", "L", "R", "F", "S", "SH", 
    "TH", "DH", "V", "Z", "ZH", "W", "Y", "HH", "CH", "JH"
};

std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> result(num);
    if (num == 1) {
        result[0] = start;
        return result;
    }
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        result[i] = start + step * i;
    }
    return result;
}

class IIRFilter {
public:
    std::vector<double> b;
    std::vector<double> a;
    std::vector<double> z;
    int order;

    IIRFilter(const std::vector<double>& b_coeffs, const std::vector<double>& a_coeffs) 
        : b(b_coeffs), a(a_coeffs) {
        order = std::max(b.size(), a.size()) - 1;
        z.assign(order, 0.0);
    }

    double process(double input) {
        double output = b[0] * input + z[0];
        for (int i = 0; i < order - 1; ++i) {
            z[i] = b[i + 1] * input - a[i + 1] * output + z[i + 1];
        }
        return output; 
    }
    
    void reset() {
        std::fill(z.begin(), z.end(), 0.0);
    }
};

std::vector<double> lfilter(const std::vector<double>& b, const std::vector<double>& a, const std::vector<double>& x) {
    int n = x.size();
    std::vector<double> y(n, 0.0);
    int nb = b.size();
    int na = a.size();
    std::vector<double> w(std::max(nb, na), 0.0);

    for (int i = 0; i < n; ++i) {
        for (int k = std::max(nb, na) - 1; k > 0; --k) {
            w[k] = w[k - 1];
        }
        w[0] = x[i];
        
        double acc = 0.0;
        for (int k = 0; k < nb; ++k) {
            acc += b[k] * w[k];
        }
        for (int k = 1; k < na; ++k) {
            acc -= a[k] * y[i - k];
        }
        y[i] = 0.0;
        for(int k=0; k<nb; ++k) {
            if(i-k >= 0) y[i] += b[k] * x[i-k];
        }
        for(int k=1; k<na; ++k) {
            if(i-k >= 0) y[i] -= a[k] * y[i-k];
        }
        if (a[0] != 1.0 && a[0] != 0.0) {
            y[i] /= a[0];
        }
    }
    return y;
}

std::vector<double> filtfilt(const std::vector<double>& b, const std::vector<double>& a, const std::vector<double>& x) {
    int padlen = 3 * (std::max(b.size(), a.size()) - 1);
    std::vector<double> x_padded = x;
    
    for(int i=0; i<padlen; ++i) {
        if(!x.empty()) {
            x_padded.insert(x_padded.begin(), x[0]);
            x_padded.push_back(x.back());
        }
    }

    std::vector<double> y_fwd = lfilter(b, a, x_padded);
    std::reverse(y_fwd.begin(), y_fwd.end());
    std::vector<double> y_bwd = lfilter(b, a, y_fwd);
    std::reverse(y_bwd.begin(), y_bwd.end());

    if (y_bwd.size() > 2 * padlen) {
        return std::vector<double>(y_bwd.begin() + padlen, y_bwd.end() - padlen);
    }
    return y_bwd;
}

void design_butterworth(double sample_rate, double freq, const std::string& type, int order, 
                        std::vector<double>& b, std::vector<double>& a) {
    
    b.clear(); a.clear();
    
    double w0 = 2.0 * M_PI * freq / sample_rate;
    double cos_w0 = std::cos(w0);
    double sin_w0 = std::sin(w0);
    
    double alpha = sin_w0 / 2.0;

    double b0, b1, b2, a0, a1, a2;

    if (type == "low") {
        b0 = (1 - cos_w0) / 2;
        b1 = 1 - cos_w0;
        b2 = (1 - cos_w0) / 2;
        a0 = 1 + alpha;
        a1 = -2 * cos_w0;
        a2 = 1 - alpha;
    } else if (type == "high") {
        b0 = (1 + cos_w0) / 2;
        b1 = -(1 + cos_w0);
        b2 = (1 + cos_w0) / 2;
        a0 = 1 + alpha;
        a1 = -2 * cos_w0;
        a2 = 1 - alpha;
    } else if (type == "band") {
        b0 = alpha;
        b1 = 0;
        b2 = -alpha;
        a0 = 1 + alpha;
        a1 = -2 * cos_w0;
        a2 = 1 - alpha;
    } else {
        b0 = (1 - cos_w0) / 2; b1 = 1 - cos_w0; b2 = (1 - cos_w0) / 2;
        a0 = 1 + alpha; a1 = -2 * cos_w0; a2 = 1 - alpha;
    }

    b = {b0/a0, b1/a0, b2/a0};
    a = {1.0, a1/a0, a2/a0};
}

struct PhonemeData {
    double f1 = 0.0, f2 = 0.0, f3 = 0.0, f4 = 0.0;
    double length = 0.14;
    bool voiced = false;
    bool is_silence = false;
};

class Voice {
public:
    std::string name;
    std::string description;
    std::map<std::string, PhonemeData> phonemes;

    Voice(const std::string& n, const std::string& desc = "") : name(n), description(desc) {}

    PhonemeData get_phoneme_data(const std::string& phoneme) {
        std::string base_ph = phoneme;
        bool is_final = false;
        if (base_ph.size() > 6 && base_ph.substr(base_ph.size() - 6) == "_FINAL") {
            base_ph = base_ph.substr(0, base_ph.size() - 6);
            is_final = true;
        }

        auto it = phonemes.find(base_ph);
        if (it == phonemes.end()) {
            auto sil_it = phonemes.find("SIL");
            if (sil_it != phonemes.end()) return sil_it->second;
            return PhonemeData();
        }

        PhonemeData data = it->second;
        if (is_final && VOWELS.find(base_ph) != VOWELS.end()) {
            data.length = std::min(data.length * 1.4, 0.35);
        }
        return data;
    }
};

class DefaultVoice : public Voice {
public:
    DefaultVoice() : Voice("Default", "Built-in robotic voice") {
        auto add = [&](const std::string& p, double f1, double f2, double f3, double f4, double len, bool voiced, bool sil=false) {
            phonemes[p] = {f1, f2, f3, f4, len, voiced, sil};
        };

        add("AH", 700, 1100, 2400, 115, 0.14, true);
        add("AE", 650, 1250, 2500, 115, 0.14, true);
        add("AA", 620, 1180, 2550, 115, 0.14, true);
        add("AO", 550, 850, 2400, 115, 0.14, true);
        add("EH", 530, 1700, 2450, 115, 0.14, true);
        add("EY", 400, 2100, 2800, 115, 0.14, true);
        add("IH", 420, 1950, 2500, 115, 0.14, true);
        add("IY", 300, 2250, 3000, 115, 0.14, true);
        add("OW", 450, 900, 2350, 115, 0.14, true);
        add("UH", 400, 650, 2400, 115, 0.14, true);
        add("UW", 330, 900, 2200, 115, 0.14, true);
        add("ER", 480, 1180, 1650, 115, 0.14, true);
        add("M",  350, 1050, 2250, 115, 0.12, true);
        add("N",  320, 1150, 2450, 115, 0.12, true);
        add("NG", 280, 950, 2350, 115, 0.12, true);
        add("L",  400, 1150, 2450, 115, 0.12, true);
        add("R",  450, 1250, 1500, 115, 0.12, true);
        add("DH", 380, 1650, 2450, 115, 0.12, true);
        add("V",  380, 1550, 2450, 115, 0.12, true);
        add("Z",  380, 1750, 2450, 115, 0.12, true);
        add("ZH", 380, 1450, 2250, 115, 0.12, true);
        add("W",  350, 700, 2450, 115, 0.12, true);
        add("Y",  350, 2050, 2650, 115, 0.12, true);
        add("JH", 400, 1650, 2450, 115, 0.12, true);
        
        add("B",  0, 0, 0, 0, 0.068, false);
        add("D",  0, 0, 0, 0, 0.068, false);
        add("G",  0, 0, 0, 0, 0.068, false);
        add("P",  0, 0, 0, 0, 0.068, false);
        add("T",  0, 0, 0, 0, 0.068, false);
        add("K",  0, 0, 0, 0, 0.068, false);
        add("F",  0, 0, 0, 0, 0.125, false);
        add("S",  0, 0, 0, 0, 0.125, false);
        add("SH", 0, 0, 0, 0, 0.125, false);
        add("TH", 0, 0, 0, 0, 0.125, false);
        add("HH", 0, 0, 0, 0, 0.125, false);
        add("CH", 0, 0, 0, 0, 0.068, false);
        
        add("SIL", 0, 0, 0, 0, 0.19, false, true);
    }
};

class VoiceRegistry {
public:
    std::map<std::string, Voice*> voices;
    Voice* current_voice;

    VoiceRegistry() {
        voices["Default"] = new DefaultVoice();
        current_voice = voices["Default"];
    }

    ~VoiceRegistry() {
        for (auto& pair : voices) {
            delete pair.second;
        }
    }

    bool set_current_voice(const std::string& name) {
        auto it = voices.find(name);
        if (it != voices.end()) {
            current_voice = it->second;
            return true;
        }
        return false;
    }
};

VoiceRegistry VOICE_REGISTRY;

std::vector<std::string> parse_phoneme_input(const std::string& text) {
    std::vector<std::string> phonemes;
    std::istringstream iss(text);
    std::string token;
    
    while (iss >> token) {
        std::transform(token.begin(), token.end(), token.begin(), ::toupper);
        if (VALID_PHONEMES.find(token) != VALID_PHONEMES.end()) {
            phonemes.push_back(token);
        } else {
            phonemes.push_back("SIL");
        }
    }

    if (phonemes.empty()) return {"SIL", "SIL"};
    if (phonemes.front() != "SIL") phonemes.insert(phonemes.begin(), "SIL");
    if (phonemes.back() != "SIL") phonemes.push_back("SIL");

    return phonemes;
}

struct PhonemeSpec {
    std::string phoneme;
    double duration;
    double overlap;
    std::vector<double> pitch_contour;
    double f1, f2, f3;
    bool voiced;
};

std::vector<PhonemeSpec> parse_phoneme_spec(const std::string& text, Voice* voice, double pitch_base = 115.0) {
    std::vector<PhonemeSpec> specs;
    std::string segment;
    std::stringstream pipe_stream(text);

    while (std::getline(pipe_stream, segment, '|')) {
        segment.erase(0, segment.find_first_not_of(" \t\r\n"));
        segment.erase(segment.find_last_not_of(" \t\r\n") + 1);
        if (segment.empty() || segment[0] == '#') continue;

        std::istringstream line_ss(segment);
        std::string ph_name;
        line_ss >> ph_name;
        std::transform(ph_name.begin(), ph_name.end(), ph_name.begin(), ::toupper);

        if (VALID_PHONEMES.find(ph_name) == VALID_PHONEMES.end()) continue;

        double duration, overlap;
        line_ss >> duration >> overlap;
        
        std::vector<double> pitch_points;
        double p;
        while (line_ss >> p) {
            pitch_points.push_back(p);
        }

        if (pitch_points.empty()) {
            pitch_points.push_back(pitch_base);
        }
        if (pitch_points.size() > 8) pitch_points.resize(8);

        duration = std::max(0.01, std::min(2.0, duration));
        overlap = std::max(0.0, std::min(0.5, overlap));

        PhonemeData ph_data = voice->get_phoneme_data(ph_name);
        
        std::set<std::string> unvoiced_set = {"SIL", "B", "D", "G", "P", "T", "K", "F", "S", "SH", "TH", "HH", "CH"};
        bool is_voiced = (unvoiced_set.find(ph_name) == unvoiced_set.end());

        specs.push_back({
            ph_name, duration, overlap, pitch_points,
            ph_data.f1, ph_data.f2, ph_data.f3, is_voiced
        });
    }
    return specs;
}

std::vector<PhonemeSpec> phonemes_to_spec(const std::vector<std::string>& phonemes, Voice* voice, double pitch_base = 115.0) {
    std::vector<PhonemeSpec> specs;
    std::set<std::string> unvoiced_set = {"SIL", "B", "D", "G", "P", "T", "K", "F", "S", "SH", "TH", "HH", "CH"};

    for (size_t i = 0; i < phonemes.size(); ++i) {
        const std::string& ph = phonemes[i];
        PhonemeData ph_data = voice->get_phoneme_data(ph);
        double duration = ph_data.length;
        double overlap = (VOWELS.find(ph) != VOWELS.end() && i < phonemes.size() - 1) ? 0.050 : 0.025;

        std::vector<double> pitch;
        if (ph == "SIL") {
            pitch = {0.0};
        } else if (VOWELS.find(ph) != VOWELS.end()) {
            if (i == phonemes.size() - 2) {
                pitch = {pitch_base * 0.95, pitch_base * 0.90};
            } else if (i == 1) {
                pitch = {pitch_base * 1.05, pitch_base * 1.10};
            } else {
                pitch = {pitch_base};
            }
        } else {
            pitch = {ph_data.voiced ? pitch_base : 0.0};
        }

        bool is_voiced = (unvoiced_set.find(ph) == unvoiced_set.end());

        specs.push_back({
            ph, duration, overlap, pitch,
            ph_data.f1, ph_data.f2, ph_data.f3, is_voiced
        });
    }
    return specs;
}

class FormantSynthesizer {
private:
    int fs;
    Voice* voice;
    double nyquist;
    int ref_fs = 44100;
    std::mt19937 rng;
    std::normal_distribution<double> dist;

public:
    FormantSynthesizer(Voice* v, int sample_rate = 44100) 
        : fs(sample_rate), voice(v), dist(0.0, 1.0) {
        nyquist = fs / 2.0;
        rng.seed(42);
    }

    double scale_freq(double freq) {
        if (freq <= 0) return 0.0;
        double scaled = freq * (static_cast<double>(fs) / ref_fs);
        return std::max(20.0, std::min(scaled, nyquist * 0.99));
    }

    std::vector<double> generate_glottal_pulse_train_contour(double duration, const std::vector<double>& pitch_contour) {
        int n_samples = static_cast<int>(duration * fs);
        std::vector<double> signal(n_samples, 0.0);
        double t = 0.0;

        std::vector<double> p_contour = pitch_contour;
        if (p_contour.empty() || std::all_of(p_contour.begin(), p_contour.end(), [](double p){ return p == 0.0; })) {
            p_contour = {115.0};
        }

        int num_points = p_contour.size();

        while (t < duration) {
            double t_norm = std::min(1.0, t / duration);
            double f0 = 0.0;

            if (num_points == 1) {
                f0 = p_contour[0];
            } else {
                double contour_pos = t_norm * (num_points - 1);
                int idx_floor = static_cast<int>(contour_pos);
                double frac = contour_pos - idx_floor;

                if (idx_floor >= num_points - 1) {
                    f0 = p_contour.back();
                } else {
                    f0 = p_contour[idx_floor] * (1 - frac) + p_contour[idx_floor + 1] * frac;
                }
            }

            double nyquist_limit = (fs / 2.0) - 500;
            f0 = std::max(50.0, std::min(nyquist_limit, f0));
            double period_samples = fs / f0;
            int pulse_len = static_cast<int>(period_samples * 0.6);
            if (pulse_len < 8) pulse_len = 8;

            std::vector<double> pulse(pulse_len, 0.0);
            int open_len = std::max(4, static_cast<int>(pulse_len * 0.4));
            
            std::vector<double> cos_table = linspace(0, M_PI, open_len);
            for(int i=0; i<open_len; ++i) {
                pulse[i] = -0.5 * (1 - std::cos(cos_table[i]));
            }

            if (pulse_len > open_len) {
                int close_len = pulse_len - open_len;
                std::vector<double> exp_table = linspace(0, 5, close_len);
                for(int i=0; i<close_len; ++i) {
                    pulse[open_len + i] = -0.1 * std::exp(-exp_table[i]);
                }
            }

            int start = static_cast<int>(t * fs);
            int end = std::min(start + pulse_len, n_samples);

            if (end > start) {
                int copy_len = std::min(pulse_len, end - start);
                for(int i=0; i<copy_len; ++i) {
                    signal[start + i] += pulse[i] * 0.6;
                }
            }

            t += period_samples / fs;
        }

        double peak = 0.0;
        for(double s : signal) peak = std::max(peak, std::abs(s));
        if (peak > 0.1) {
            double scale = 0.6 / peak;
            for(double& s : signal) s *= scale;
        }

        return signal;
    }

    std::vector<double> generate_shaped_noise(double duration, const std::string& phoneme, double intensity = 0.25, double f0 = 115.0) {
        int n_samples = static_cast<int>(duration * fs);
        std::vector<double> noise(n_samples);
        for(int i=0; i<n_samples; ++i) noise[i] = dist(rng);

        std::vector<double> b, a;

        if (phoneme == "S") {
            double low_cut  = scale_freq(6500);   
            double high_cut = scale_freq(12000);  

            design_butterworth(fs, low_cut, "high", 4, b, a);
            noise = lfilter(b, a, noise);

            design_butterworth(fs, high_cut, "low", 4, b, a);
            noise = lfilter(b, a, noise);

            for (double& n : noise) n *= 1.6;
        }
        else if (phoneme == "SH" || phoneme == "ZH") {

            double low = scale_freq(2500);
            double high = scale_freq(6000);

            design_butterworth(fs, low, "high", 2, b, a);
            noise = lfilter(b, a, noise);

            design_butterworth(fs, high, "low", 2, b, a);
            noise = lfilter(b, a, noise);

            for(int i=0;i<n_samples;i++)
                noise[i] *= 0.07;   
        }
        else if (phoneme == "TH") {

            double low = scale_freq(800);
            double high = scale_freq(4500);

            design_butterworth(fs, low, "high", 2, b, a);
            noise = lfilter(b, a, noise);

            design_butterworth(fs, high, "low", 2, b, a);
            noise = lfilter(b, a, noise);

            for(int i=0;i<n_samples;i++)
                noise[i] *= 0.03;   
        }
        else if (phoneme == "F") {

            double low = scale_freq(1000);
            double high = scale_freq(7000);

            design_butterworth(fs, low, "high", 2, b, a);
            noise = lfilter(b, a, noise);

            design_butterworth(fs, high, "low", 2, b, a);
            noise = lfilter(b, a, noise);

            for(int i=0;i<n_samples;i++)
                noise[i] *= 0.05;   
        }
        else if (phoneme == "HH") {
            double cutoff_f = scale_freq(1800);
            design_butterworth(fs, cutoff_f, "low", 2, b, a);
            noise = lfilter(b, a, noise);
            for(int i=0; i<n_samples; ++i) noise[i] += dist(rng) * 0.9;
        }
        else if (phoneme == "Z") {

            double low = scale_freq(3500);
            double high = scale_freq(8000);

            design_butterworth(fs, low, "high", 2, b, a);
            noise = lfilter(b, a, noise);

            design_butterworth(fs, high, "low", 2, b, a);
            noise = lfilter(b, a, noise);

            for(int i=0; i<n_samples; ++i) {
                double phase = 2 * M_PI * f0 * i / fs;
                double voicing = (std::sin(phase)
                                 + 0.5 * std::sin(2 * phase)
                                 + 0.25 * std::sin(3 * phase)) * 0.08;

                noise[i] = noise[i] * 0.7 + voicing * 0.3;
            }
        }
        else if (phoneme == "V") {

            double low = scale_freq(1000);
            double high = scale_freq(7000);

            design_butterworth(fs, low, "high", 2, b, a);
            noise = lfilter(b, a, noise);

            design_butterworth(fs, high, "low", 2, b, a);
            noise = lfilter(b, a, noise);

            for(int i=0; i<n_samples; ++i) {
                double phase = 2 * M_PI * f0 * i / fs;
                double voicing = (std::sin(phase)
                                 + 0.5 * std::sin(2 * phase)
                                 + 0.25 * std::sin(3 * phase)) * 0.12;

                noise[i] = noise[i] * 0.5 + voicing * 0.5;
            }
        }
        else if (phoneme == "DH") {

            double low = scale_freq(800);
            double high = scale_freq(4500);

            design_butterworth(fs, low, "high", 2, b, a);
            noise = lfilter(b, a, noise);

            design_butterworth(fs, high, "low", 2, b, a);
            noise = lfilter(b, a, noise);

            for(int i=0; i<n_samples; ++i) {
                double phase = 2 * M_PI * f0 * i / fs;
                double voicing = (std::sin(phase)
                                 + 0.5 * std::sin(2 * phase)
                                 + 0.25 * std::sin(3 * phase)) * 0.15;

                noise[i] = noise[i] * 0.4 + voicing * 0.6;
            }
        }
        else {
            double cutoff_f = scale_freq(7500);
            design_butterworth(fs, cutoff_f, "low", 2, b, a);
            noise = lfilter(b, a, noise);
        }

        double peak = 0.0;
        for(double n : noise) peak = std::max(peak, std::abs(n));
        
        if (peak < 1e-6) {
            for(int i=0; i<n_samples; ++i) noise[i] = dist(rng) * intensity * 0.7;
            peak = 1.0;
        }

        double scale = intensity / peak;
        for(double& n : noise) n *= scale;

        if (noise.size() > n_samples) noise.resize(n_samples);
        return noise;
    }

    std::vector<double> apply_formants_safe(const std::vector<double>& signal, double f1, double f2, double f3) {
        std::vector<double> out = signal;
        double bw1 = 100, bw2 = 150, bw3 = 250; 
        std::vector<double> b, a;

        auto apply_res = [&](double freq, double bw) {
            if (freq <= 50) return;
            double w0 = 2.0 * M_PI * freq / fs;
            double sin_w0 = std::sin(w0);
            double cos_w0 = std::cos(w0);
            
            double Q = freq / bw;
            double alpha = sin_w0 / (2.0 * Q);
            
            double b0 =   sin_w0 / 2.0;
            double b1 =   0.0;
            double b2 =  -sin_w0 / 2.0;
            double a0 =   1.0 + alpha;
            double a1 =  -2.0 * cos_w0;
            double a2 =   1.0 - alpha;
            
            b = {b0/a0, b1/a0, b2/a0};
            a = {1.0, a1/a0, a2/a0};
            
            out = lfilter(b, a, out);
        };

        if (f1 > 50) apply_res(f1, bw1);
        if (f2 > 50) apply_res(f2, bw2);
        if (f3 > 50) apply_res(f3, bw3);

        double peak = 0.0;
        for(double s : out) peak = std::max(peak, std::abs(s));
        if (peak > 1e-6 && peak > 3.0) {
            double scale = 2.5 / peak;
            for(double& s : out) s *= scale;
        }

        return out;
    }

    std::vector<double> synthesize_phoneme_direct(const PhonemeSpec& spec) {
        std::string ph = spec.phoneme;
        double dur = spec.duration;
        double f1 = spec.f1, f2 = spec.f2, f3 = spec.f3;
        bool voiced = spec.voiced;

        if (ph == "SIL") {
            return std::vector<double>(static_cast<int>(dur * fs), 0.0);
        }

        if (STOPS.find(ph) != STOPS.end()) {
            int n_samples = static_cast<int>(dur * fs);
            std::vector<double> out(n_samples, 0.0);

            int closure_end = static_cast<int>(n_samples * 0.80);
            for (int i = 0; i < closure_end; ++i) {
                out[i] = dist(rng) * 0.002;
            }
            int burst_start = closure_end;
            int burst_len = std::min(200, n_samples - burst_start);

            if (burst_len > 30) {
                std::vector<double> burst_noise(burst_len);
                for(int i=0; i<burst_len; ++i)
                    burst_noise[i] = dist(rng);

                double low = 1000, high = 4000;

                if (ph == "P" || ph == "B") {
                    low = 500;  high = 1500;
                }
                else if (ph == "T" || ph == "D") {
                    low = 3000; high = 5000;
                }
                else if (ph == "K" || ph == "G") {
                    low = 1500; high = 3000;
                }
                else if (ph == "CH") {
                    low = 3000; high = 6000;
                }

                std::vector<double> b,a;

                design_butterworth(fs, scale_freq(low), "high", 2, b, a);
                burst_noise = lfilter(b, a, burst_noise);

                design_butterworth(fs, scale_freq(high), "low", 2, b, a);
                burst_noise = lfilter(b, a, burst_noise);

                for(int i=0;i<burst_len;i++) {
                    double decay = std::exp(-6.0 * i / burst_len);
                    burst_noise[i] *= decay * 0.5;
                }

                for(int i=0;i<burst_len;i++)
                    out[burst_start + i] += burst_noise[i];
            }
                        
            bool is_voiced_stop = (ph == "B" || ph == "D" || ph == "G");

            if (is_voiced_stop) {

                int prevoice_len = std::min(closure_end, n_samples / 3);

                for(int i = 0; i < prevoice_len; ++i) {
                    double phase = 2 * M_PI * 100.0 * i / fs;
                    double voicing = std::sin(phase) * 0.25;   
                    out[i] += voicing;
                }

                int voice_start = burst_start + burst_len;
                int ramp_len = std::min(400, n_samples - voice_start);

                for(int i = 0; i < ramp_len; ++i) {
                    double phase = 2 * M_PI * 120.0 * i / fs;
                    double glottal = std::sin(phase) * 0.35;  
                    double fade = i / static_cast<double>(ramp_len);
                    out[voice_start + i] += glottal * fade;
                }
            }
            
            bool is_voiceless = (ph == "P" || ph == "T" || ph == "K" || ph == "CH");

            if (is_voiceless) {
                int asp_len = std::min(400, n_samples - (burst_start + burst_len));
                for(int i=0;i<asp_len;i++) {
                    double noise = dist(rng) * 0.2;
                    double decay = std::exp(-3.0 * i / asp_len);
                    out[burst_start + burst_len + i] += noise * decay;
                }
            }

            if (!is_voiced_stop)
                for(double& o : out) o *= 0.85;
            return out;
        }

        std::vector<double> output;
        if (!voiced) {
            double intensity = 0.32;
            if (ph == "S") intensity = 1.28;
            else if (ph == "SH") intensity = 0.64;
            else if (ph == "F" || ph == "TH") intensity = 0.32;

            std::vector<double> source = generate_shaped_noise(dur, ph, intensity, spec.pitch_contour.empty() ? 115.0 : spec.pitch_contour[0]);
            if (f1 > 50) {
                if (ph == "S" || ph == "SH") {
                    source = apply_formants_safe(source, f1*0.7, f2*0.7, f3*0.7);
                } else {
                    source = apply_formants_safe(source, f1, f2, f3);
                }
            } else {
                for(double& s : source) s *= 0.45;
            }
            output = source;
        } else {
            std::vector<double> source = generate_glottal_pulse_train_contour(dur, spec.pitch_contour);
            if (f1 > 50) {
                output = apply_formants_safe(source, f1, f2, f3);
            } else {
                output = source;
                for(double& o : output) o *= 0.45;
            }

            int n = output.size();
            std::vector<double> env(n, 1.0);
            double att = std::min(0.015, dur * 0.25); 
            double rel = std::min(0.030, dur * 0.40);
            int att_s = static_cast<int>(att * fs);
            int rel_s = static_cast<int>(rel * fs);

            if (att_s > 0) {
                std::vector<double> ramp = linspace(0, 1, att_s);
                for(int i=0; i<att_s && i<n; ++i) env[i] = ramp[i];
            }
            if (rel_s > 0) {
                std::vector<double> ramp = linspace(1, 0.05, rel_s);
                for(int i=0; i<rel_s && (n-1-i)>=0; ++i) env[n-1-i] = ramp[i];
            }

            for(int i=0; i<n; ++i) output[i] *= env[i];
            for(double& o : output) o = std::tanh(o * 1.15) * 0.93;
        }

        for(double& o : output) o *= 0.82;
        return output;
    }

    std::vector<double> synthesize_from_specs(const std::vector<PhonemeSpec>& specs) {
        if (specs.empty()) return {};

        double total_duration = 0.0;
        for (const auto& spec : specs) total_duration += spec.duration;
        for (size_t i = 0; i < specs.size() - 1; ++i) {
            total_duration -= std::min(specs[i].overlap, specs[i].duration);
        }

        int total_samples = static_cast<int>(total_duration * fs) + 10;
        std::vector<double> output(total_samples, 0.0);
        int current_pos = 0;

        for (size_t i = 0; i < specs.size(); ++i) {
            std::vector<double> phoneme_audio = synthesize_phoneme_direct(specs[i]);
            int phoneme_samples = phoneme_audio.size();
            double overlap_dur = specs[i].overlap;
            if (i == specs.size() - 1) overlap_dur = 0.0;

            int overlap_samples = std::min({
                static_cast<int>(overlap_dur * fs),
                phoneme_samples - 1,
                static_cast<int>(specs[i].duration * fs * 0.5)
            });

            int end_pos = current_pos + phoneme_samples;
            if (end_pos > output.size()) {
                output.resize(end_pos + 1000, 0.0);
            }

            for(int k=0; k<phoneme_samples; ++k) {
                if(current_pos + k < output.size()) {
                    output[current_pos + k] += phoneme_audio[k];
                }
            }
            current_pos += (phoneme_samples - overlap_samples);
        }

        int actual_length = std::min(current_pos, static_cast<int>(output.size()));
        std::vector<double> audio(output.begin(), output.begin() + actual_length);

        for(double& a : audio) a = std::tanh(a * 1.25) * 0.94;

        double cutoff_f = scale_freq(5000);
        std::vector<double> b, a;
        design_butterworth(fs, cutoff_f, "low", 2, b, a);
        audio = lfilter(b, a, audio);

        return audio;
    }
};

std::vector<uint8_t> generate_wav_bytes(const std::vector<double>& audio, int sr) {
    std::vector<int16_t> audio_int16(audio.size());
    for (size_t i = 0; i < audio.size(); ++i) {
        double val = audio[i] * 32767.0;
        if (val > 32767.0) val = 32767.0;
        if (val < -32768.0) val = -32768.0;
        audio_int16[i] = static_cast<int16_t>(val);
    }

    std::vector<uint8_t> buffer;
    int num_channels = 1;
    int bits_per_sample = 16;
    int byte_rate = sr * num_channels * bits_per_sample / 8;
    int block_align = num_channels * bits_per_sample / 8;
    int data_size = audio_int16.size() * 2;

    auto write_u32 = [&](uint32_t val) {
        buffer.push_back(val & 0xFF);
        buffer.push_back((val >> 8) & 0xFF);
        buffer.push_back((val >> 16) & 0xFF);
        buffer.push_back((val >> 24) & 0xFF);
    };
    auto write_u16 = [&](uint16_t val) {
        buffer.push_back(val & 0xFF);
        buffer.push_back((val >> 8) & 0xFF);
    };

    buffer.insert(buffer.end(), {'R', 'I', 'F', 'F'});
    write_u32(36 + data_size);
    buffer.insert(buffer.end(), {'W', 'A', 'V', 'E'});

    buffer.insert(buffer.end(), {'f', 'm', 't', ' '});
    write_u32(16);
    write_u16(1);
    write_u16(num_channels);
    write_u32(sr);
    write_u32(byte_rate);
    write_u16(block_align);
    write_u16(bits_per_sample);

    buffer.insert(buffer.end(), {'d', 'a', 't', 'a'});
    write_u32(data_size);

    for (int16_t sample : audio_int16) {
        buffer.push_back(sample & 0xFF);
        buffer.push_back((sample >> 8) & 0xFF);
    }

    return buffer;
}

std::vector<uint8_t> synthesize_phoneme_mode(const std::string& text, int sample_rate = 44100, Voice* voice = nullptr) {
    if (!voice) voice = VOICE_REGISTRY.current_voice;
    std::vector<std::string> phonemes = parse_phoneme_input(text);
    std::vector<PhonemeSpec> specs = phonemes_to_spec(phonemes, voice);
    FormantSynthesizer synth(voice, sample_rate);
    std::vector<double> audio = synth.synthesize_from_specs(specs);
    return generate_wav_bytes(audio, sample_rate);
}

std::vector<uint8_t> synthesize_spec_mode(const std::string& text, int sample_rate = 44100, Voice* voice = nullptr) {
    if (!voice) voice = VOICE_REGISTRY.current_voice;
    std::vector<PhonemeSpec> specs = parse_phoneme_spec(text, voice);
    FormantSynthesizer synth(voice, sample_rate);
    std::vector<double> audio = synth.synthesize_from_specs(specs);
    return generate_wav_bytes(audio, sample_rate);
}

}

int main(int argc, char* argv[]) {
    using namespace DaSpeakSynth;
    
    std::string spec_input;
    std::string phon_input;
    double volume_db = 0.0;
    double pitch_base = 115.0;  
    int sample_rate = 44100;
    std::string output_file = "output.wav";
    std::string voice_name = "Default";
    bool show_help = false;
    
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            show_help = true;
            continue;
        }
        
        if (arg.rfind("-", 0) == 0 && arg.find('=') != std::string::npos) {
            size_t eq_pos = arg.find('=');
            std::string key = arg.substr(0, eq_pos);
            std::string value = arg.substr(eq_pos + 1);
            
            if (value.size() >= 2 && 
                ((value.front() == '"' && value.back() == '"') || 
                 (value.front() == '\'' && value.back() == '\''))) {
                value = value.substr(1, value.size() - 2);
            }
            
            if (key == "-spec") {
                spec_input = value;
            } else if (key == "-phon") {
                phon_input = value;
            } else if (key == "-v") {
                try {
                    volume_db = std::stod(value);
                } catch (...) {
                    std::cerr << "Error: Invalid volume value '" << value << "'\n";
                    return 1;
                }
            } else if (key == "-p") {  
                try {
                    pitch_base = std::stod(value);
                    if (pitch_base < 50.0 || pitch_base > 500.0) {
                        std::cerr << "Warning: Pitch " << pitch_base 
                                  << " Hz outside recommended range (50-500 Hz)\n";
                    }
                } catch (...) {
                    std::cerr << "Error: Invalid pitch value '" << value << "'\n";
                    return 1;
                }
            } else if (key == "-r") {
                try {
                    sample_rate = std::stoi(value);
                    if (sample_rate < 8000 || sample_rate > 192000) {
                        std::cerr << "Warning: Sample rate " << sample_rate 
                                  << " Hz is outside typical range (8000-192000)\n";
                    }
                } catch (...) {
                    std::cerr << "Error: Invalid sample rate '" << value << "'\n";
                    return 1;
                }
            } else if (key == "-o") {
                output_file = value;
                if (output_file.size() < 4 || 
                    output_file.substr(output_file.size() - 4) != ".wav") {
                    output_file += ".wav";
                }
            } else if (key == "-voice") {
                voice_name = value;
            } else {
                std::cerr << "Warning: Unknown option '" << key << "'\n";
            }
        } else if (arg.rfind("-", 0) == 0) {
            std::cerr << "Warning: Malformed option '" << arg << "' (expected key=value)\n";
        } else {
            std::cerr << "Warning: Unexpected argument '" << arg << "'\n";
        }
    }
    
    if (show_help || (spec_input.empty() && phon_input.empty())) {
        std::cout << "DaSpeak Formant Synthesizer v1.0\n"
                  << "Usage: " << argv[0] << " [options]\n\n"
                  << "Options:\n"
                  << "  -spec=\"<text>\"   Phoneme specifications (pipe-separated)\n"
                  << "                     Format: PHONEME DURATION OVERLAP PITCH...\n"
                  << "                     Example: \"HH 0.12 0.015 95|EH 0.14 0.018 105 110\"\n"
                  << "  -phon=\"<text>\"   Space-separated phoneme names\n"
                  << "                     Example: \"HH EH L OW\"\n"
                  << "  -p=<Hz>          Base pitch frequency in Hz (default: 115)\n"
                  << "                     Used for -phon mode; fallback for -spec if pitch omitted\n"
                  << "  -v=<dB>          Volume adjustment in decibels (default: 0)\n"
                  << "  -r=<rate>        Sample rate in Hz (default: 44100)\n"
                  << "  -o=<file>        Output WAV filename (default: output.wav)\n"
                  << "  -voice=<name>    Voice preset name (default: Default)\n"
                  << "  -h, --help       Show this help message\n\n"
                  << "Note: Either -spec or -phon must be provided.\n"
                  << "      -spec allows fine control over duration, overlap, and pitch contour.\n"
                  << "      -phon uses default timing and pitch for each phoneme.\n";
        return (show_help ? 0 : 1);
    }
    
    if (!spec_input.empty() && !phon_input.empty()) {
        std::cerr << "Error: Cannot specify both -spec and -phon. Choose one.\n";
        return 1;
    }
    if (spec_input.empty() && phon_input.empty()) {
        std::cerr << "Error: Must specify either -spec or -phon.\n";
        return 1;
    }
    
    if (!VOICE_REGISTRY.set_current_voice(voice_name)) {
        std::cerr << "Warning: Voice '" << voice_name << "' not found, using Default\n";
        VOICE_REGISTRY.set_current_voice("Default");
    }
    Voice* voice = VOICE_REGISTRY.current_voice;
    
    std::vector<double> audio;
    try {
        if (!spec_input.empty()) {
            std::vector<PhonemeSpec> specs = parse_phoneme_spec(spec_input, voice, pitch_base);
            if (specs.empty()) {
                std::cerr << "Error: No valid phoneme specifications parsed\n";
                return 1;
            }
            FormantSynthesizer synth(voice, sample_rate);
            audio = synth.synthesize_from_specs(specs);
        } else {
            FormantSynthesizer synth(voice, sample_rate);
            std::vector<std::string> phonemes = parse_phoneme_input(phon_input);
            std::vector<PhonemeSpec> specs = phonemes_to_spec(phonemes, voice, pitch_base);
            audio = synth.synthesize_from_specs(specs);
        }
    } catch (const std::exception& e) {
        std::cerr << "Synthesis error: " << e.what() << "\n";
        return 1;
    }
    
    if (audio.empty()) {
        std::cerr << "Error: Synthesis produced no audio data\n";
        return 1;
    }
    
    if (volume_db != 0.0) {
        double gain = std::pow(10.0, volume_db / 20.0);
        for (double& sample : audio) {
            sample *= gain;
            sample = std::tanh(sample);
        }
    }
    
    std::vector<uint8_t> wav_data = generate_wav_bytes(audio, sample_rate);
    
    std::ofstream outfile(output_file, std::ios::binary);
    if (!outfile) {
        std::cerr << "Error: Cannot open output file '" << output_file << "'\n";
        return 1;
    }
    outfile.write(reinterpret_cast<const char*>(wav_data.data()), wav_data.size());
    outfile.close();
    
    std::cout << "Successfully synthesized " << audio.size() << " samples (" 
              << std::fixed << std::setprecision(2) 
              << static_cast<double>(audio.size()) / sample_rate << " seconds)\n"
              << "Output written to: " << output_file << "\n";
    
    return 0;
}
