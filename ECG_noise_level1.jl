
"""
ECG Noise

Задача: удаление помех и шумов (подготовка к выделению QRS-комплексов)

1. Применить фильтры (для каждого файла можно выбрать один или несколько каналов):
    1.1. Выбрать файл 50hz1 - применить режекторный фильтр 50 Гц (либо заграждающий в полосе 40-60 Гц).
    1.2. Выбрать файл lf* - применить ВЧ-фильтр изолинии, чтобы убрать дрейф (отфильтровать частоты ниже 0.1 Гц).
    1.3. Выбрать файл hf* - применить НЧ-фильтр, чтобы убрать миографическую помеху. (Либо фильтр медианы в скользящем окне)
    1.4. Выбрать файл step* - выделить непригодные для анализа участки с резкими скачками 
    (применить последовательно предыдущие фильтры + компаратор, порог подобрать вручную)
    1.5. Выбрать файл chatter* - выделить непригодные для анализа участки с резкими скачками 
    (применить последовательно предыдущие фильтры + компаратор, порог подобрать вручную)

2. Для каждого файла выше - отобразить на графике участок сигнала с фильтрацией и без. 
Заголовок окна графиков - имя записи и название отведения.
Окно с графиками сохранить в виде изображений.

3. Сделать выводы об эффективности фильтрации (подавление шумов и сохранение формы полезного сигнала).
Какую комбинацию фильтров нужно использовать для фильтрации всех подобных помех?
"""
#=
   - Прочитать фрагмент сигнала и вывести на график
   - Отфильтровать предложенными способами
   - Подобрать параметры фильтров так, чтобы решить исходную задачу=#
using Plots, Dates, CSV, DataFrames, FFTW, DSP
include("readfiles.jl")
istart = 2000
len =  4000
samples = istart : istart + len - 1
channels, fs, timestart, units = readbin("data/ECG Noise/50hz1.bin", samples)#прочитать небольшой фрагмент зашумленого сигнала
channels = DataFrame(channels)
ecg1 = channels[:, :C2F] 
plot_raw = plot(ecg1, title = "50hz1.bin_C2F")

function create_filter(ftype::FilterType, coefs::FilterCoefficients)
    df = digitalfilter(ftype, coefs) |> DF2TFilter
    a = coefa(df.coef)
    b = coefb(df.coef)
    return a, b
end
#        rejection filter 
# #Fc = [40,60] # 50hz
fcut_bs1 = 40
fcut_bs2 = 60
a_bs, b_bs = create_filter(Bandstop(40, 60; fs), Butterworth(2))
ecgf = copy(ecg1)
label = ""
DSP.filt!(ecgf, b_bs, a_bs, ecgf)
label *= "+ bandstop $fcut_bs1-$fcut_bs2 Hz"
ecgf[1:50] .= NaN
plot_filter = plot(ecgf)
p_rejection=plot(plot_raw, plot_filter, layout = (2, 1), label = ["raw" label])
savefig(p_rejection,"ECG_1.1.png")

#ВЧ-фильтр -> убрать дрейф изолинии
channels, fs, timestart, units = readbin("data/ECG Noise/lf2.bin", samples) 
channels = DataFrame(channels)
ecg2 = channels[:, :C2L] 
plot_raw=plot(ecg2, title = "lf2.bin_C2L")

fcut_hp = 1
a_hp, b_hp = create_filter(Highpass(2*fcut_hp/fs), Butterworth(2))
ecgf = copy(ecg2)
DSP.filt!(ecgf, b_hp, a_hp, ecgf)
label = "+ highpass $fcut_hp Hz"
plot_filter=plot(ecgf)
p_highpass=plot(plot_raw, plot_filter,  layout = (2, 1), label = ["raw" label])
savefig(p_highpass,"ECG_1.2.png")

#НЧ-фильтр, чтобы убрать миографическую помеху
channels, fs, timestart, units = readbin("data/ECG Noise/hf3.bin", samples) 
channels = DataFrame(channels)
ecg3 = channels[:, :C1] 
plot_raw=plot(ecg3, title = "hf3.bin_C1")

# Fc = 35 # 35hz
fcut_lp = 5
a_lp, b_lp = create_filter(Lowpass(2*fcut_lp/fs), Butterworth(2))
ecgf = copy(ecg3)
DSP.filt!(ecgf, b_lp, a_lp, ecgf)
label = "+ lowpass $fcut_lp Hz"
plot_filter=plot(ecgf)
p_lowpass=plot(plot_raw, plot_filter,  layout = (2, 1), label = ["raw" label])
savefig(p_lowpass,"ECG_1.3.png")

#1.4. Выбрать файл step* - выделить непригодные для анализа участки с резкими скачками 
#   (применить последовательно предыдущие фильтры + компаратор, порог подобрать вручную)
istart = 1000
len =  10000
samples = istart : istart + len - 1
channels, fs, timestart, units = readbin("data/ECG Noise/step5.bin", samples)#прочитать небольшой фрагмент зашумленого сигнала
channels = DataFrame(channels)
ecg4 = channels[:, :C1] 
plot_raw = plot(ecg4, title = "step5.bin_C1")

#rejection filter 
fcut_bs1 = 40
fcut_bs2 = 60
ecg_r = filtfilt(digitalfilter(Bandstop(fcut_bs1, fcut_bs2; fs), Butterworth(2)), ecg4)
label_r = "+ Bandstop $fcut_bs1-$fcut_bs2 Hz"
plot_r = plot(ecg_r)

# ВЧ-фильтр (последовательно)
fcut_hp = 1
ecg_hp = filtfilt(digitalfilter(Highpass(fcut_hp; fs), Butterworth(2)), ecg_r)
label_hp = "+ highpass $fcut_hp Hz"
plot_f_hp = plot(ecg_hp)

#НЧ - фильтр (последовательно)
fcut_lp = 20
ecg_lp = filtfilt(digitalfilter(Lowpass(fcut_lp; fs), Butterworth(2)), ecg_hp)
label_lp = "+ lowpass $fcut_lp Hz"
plot_f_lp=plot(ecg_lp)

#компаратор
function detect_threshold(x)
    return x < -225 || x> 225 #ип результата отличается - Bool!
end
ecg_result = detect_threshold.(ecg_lp)

plot_result=plot(ecg_result)

#plot(plot_raw, plot_r, plot_f_hp, plot_f_lp, layout = (5, 1), label = ["raw" label_r label_hp label_lp])
result = plot(plot_raw, plot_f_lp, plot_result, layout = (3, 1), label = ["raw" "rej+hp+lp" "komporator"])
savefig(result,"ECG_1.4.png")


#1.5. Выбрать файл chatter** - выделить непригодные для анализа участки с резкими скачками 
#   (применить последовательно предыдущие фильтры + компаратор, порог подобрать вручную)
istart = 7000
len =  7000
samples = istart : istart + len - 1
channels, fs, timestart, units = readbin("data/ECG Noise/chatter4.bin", samples)#прочитать небольшой фрагмент зашумленого сигнала
channels = DataFrame(channels)
ecg5 = channels[:, :L] 
plot_raw = plot(ecg5, title = "chatter4.bin_")

#rejection filter 
fcut_bs1 = 40
fcut_bs2 = 60
ecg_r = filtfilt(digitalfilter(Bandstop(fcut_bs1, fcut_bs2; fs), Butterworth(2)), ecg5)
label_r = "+ Bandstop $fcut_bs1-$fcut_bs2 Hz"
plot_r = plot(ecg_r)

# ВЧ-фильтр (последовательно)
fcut_hp = 1.3
ecg_hp = filtfilt(digitalfilter(Highpass(fcut_hp; fs), Butterworth(2)), ecg_r)
label_hp = "+ highpass $fcut_hp Hz"
plot_f_hp = plot(ecg_hp)

#НЧ - фильтр (последовательно)
fcut_lp = 30
ecg_lp = filtfilt(digitalfilter(Lowpass(fcut_lp; fs), Butterworth(2)), ecg_hp)
label_lp = "+ lowpass $fcut_lp Hz"
plot_f_lp=plot(ecg_lp)

#компаратор
function detect_threshold(x)
    return x < -500 || x > 1000 #ип результата отличается - Bool!
end
ecg_result = detect_threshold.(ecg_lp)

plot_result=plot(ecg_result)


#plot(plot_raw, plot_r, plot_f_hp, plot_f_lp, layout = (5, 1), label = ["raw" label_r label_hp label_lp])
p_result=plot(plot_raw, plot_f_lp, plot_result, layout = (3, 1), label = ["raw" "rej+hp+lp" "komporator"])
savefig(p_result,"ECG_1.5.png")
 
