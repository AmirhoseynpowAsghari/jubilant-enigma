# -*- coding: utf-8 -*-
"""
build_pdf_fa.py -- نسخهٔ فارسی الحاقیهٔ اصلاحی arXiv:2607.02812
چینش راست‌به‌چپ با قلم B Nazanin؛ معادله‌ها به‌صورت تصویر ریاضی حروف‌چینی‌شده
(سبک LaTeX/Word، با موتور mathtext) در جریان متن قرار می‌گیرند تا ترتیب
متن فارسی/انگلیسی هرگز به‌هم نریزد.
"""
import os
import re
import pymupdf
import arabic_reshaper as ar
from bidi.algorithm import get_display
from eqrender import eq_img

HERE = os.path.dirname(os.path.abspath(__file__))
FA_TTF = os.path.join(HERE, "assets", "BNazanin.ttf")
ETA1 = os.path.join(HERE, "vendor", "eta1") + os.sep
os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
EN_TTF = "/usr/local/lib/python3.13/site-packages/matplotlib/mpl-data/fonts/ttf/DejaVuSerif.ttf"

W, H = 595.0, 842.0
MR = ML = 60.0
USABLE = W - MR - ML

faF = pymupdf.Font(fontfile=FA_TTF)
enF = pymupdf.Font(fontfile=EN_TTF)

doc = pymupdf.open()
page = doc.new_page(width=W, height=H)
Y = [70.0]
EQNO = [0]
FAD = "۰۱۲۳۴۵۶۷۸۹"


def fad(n):
    return "".join(FAD[int(d)] for d in str(n))


def fa_vis(t):
    return get_display(ar.reshape(t))


def w_fa(t, s):
    return faF.text_length(fa_vis(t), s)


def w_en(t, s):
    return enF.text_length(t, s)


def newpage():
    global page
    page = doc.new_page(width=W, height=H)
    Y[0] = 70.0


def need(h):
    if Y[0] + h > H - 60:
        newpage()


def run_metrics(r, size):
    """(width, height) of a run: ('fa',t) | ('en',t) | ('eq',latex)"""
    if r[0] == "fa":
        return w_fa(r[1], size), size
    if r[0] == "en":
        return w_en(r[1], size * 0.95), size
    _, w, h = eq_img(r[1], size)
    return w, h


def runs_width(runs, size):
    sp = faF.text_length(" ", size)
    return sum(run_metrics(r, size)[0] for r in runs) + sp * max(len(runs) - 1, 0)


def draw_runs(runs, y, size, x_start=None):
    cursor = (W - MR) if x_start is None else x_start
    sp = faF.text_length(" ", size)
    yc = y - size * 0.35
    for lang, txt in runs:
        w, h = run_metrics((lang, txt), size)
        x = cursor - w
        if lang == "eq":
            page.insert_image(pymupdf.Rect(x, yc - h / 2, x + w, yc + h / 2),
                              filename=eq_img(txt, size)[0])
        else:
            page.insert_text((x, y), fa_vis(txt) if lang == "fa" else txt,
                             fontname="fa" if lang == "fa" else "en",
                             fontfile=FA_TTF if lang == "fa" else EN_TTF,
                             fontsize=size if lang == "fa" else size * 0.95,
                             color=(0, 0, 0))
        cursor = x - sp
    return cursor


LAT = r"[A-Za-z][A-Za-z0-9_./\-]*"
ARAB = re.compile(r"[\u0600-\u06FF]")


def parse(text):
    """@@...@@ = فرمول (تصویر ریاضی)، [[...]] = متن انگلیسی اتمی، بقیه فارسی"""
    out = []
    for part in re.split(r"(@@.*?@@|\[\[.*?\]\])", text):
        if not part:
            continue
        if part.startswith("@@"):
            out.append(["eq", part[2:-2]])
            continue
        if part.startswith("[["):
            out.append(["en", part[2:-2]])
            continue
        for sub in re.split("(" + LAT + ")", part):
            if not sub:
                continue
            out.append(["en", sub] if re.fullmatch(LAT, sub) else ["fa", sub])
    for r in out:
        if r[0] == "fa" and not ARAB.search(r[1]):
            r[0] = "en"
    merged = []
    for lang, txt in out:
        if merged and merged[-1][0] == lang and lang == "fa":
            merged[-1][1] += txt
        else:
            merged.append([lang, txt])
    return merged


def wrap(runs, size):
    """خط‌بندی حریصانه؛ فرمول‌ها و ران‌های انگلیسی اتمی‌اند؛ ارتفاع هر خط هم برگردانده می‌شود"""
    sp = faF.text_length(" ", size)
    lines, cur, pend, plang = [], [], "", None
    for lang, txt in runs:
        words = [txt] if lang != "fa" else txt.split(" ")
        for word in words:
            cand = (pend + " " + word) if pend else word
            test = runs_width(cur + [(lang, cand)], size)
            if test <= USABLE or not (cur or pend):
                pend, plang = cand, lang
            else:
                if pend:
                    cur.append((plang, pend))
                if cur:
                    lines.append(cur)
                cur, pend, plang = [], word, lang
        if pend:
            cur.append((plang, pend))
            pend, plang = "", None
    if cur:
        lines.append(cur)
    out = []
    for ln in lines:
        h = max(run_metrics(r, size)[1] for r in ln)
        out.append((ln, h))
    return out


def paragraph(text, size=10.5, lead=1.75, first_gap=0.0):
    Y[0] += first_gap
    for ln, h in wrap(parse(text), size):
        lh = max(lead * size, h + 4)
        need(lh)
        draw_runs(ln, Y[0] + size + (h - size) / 2, size)
        Y[0] += lh
    Y[0] += 0.25 * size


def heading(text, size=13.5):
    need(3 * size)
    Y[0] += 0.9 * size
    draw_runs([("fa", text)], Y[0] + size, size)
    Y[0] += 1.7 * size
    page.draw_line((W - MR, Y[0] - size * 0.55), (ML, Y[0] - size * 0.55),
                   color=(0.45, 0.0, 0.0), width=0.8)
    Y[0] += 0.7 * size


def eq_display(latex, size=11.5, numbered=True):
    """معادلهٔ نمایشی وسط‌چین با شمارهٔ فارسی در حاشیهٔ چپ"""
    path, w, h = eq_img(latex, size)
    if numbered:
        EQNO[0] += 1
    need(h + 18)
    Y[0] += 6
    yc = Y[0] + h / 2
    page.insert_image(pymupdf.Rect((W - w) / 2, yc - h / 2, (W + w) / 2, yc + h / 2),
                      filename=path)
    if numbered:
        n = f"({fad(EQNO[0])})"
        page.insert_text((ML + 4, yc + 3.5), fa_vis(n), fontname="fa", fontfile=FA_TTF,
                         fontsize=10, color=(0, 0, 0))
    Y[0] += h + 10


def cell_runs(cell):
    return cell if isinstance(cell, list) else [cell]


def table(headers, rows, widths, title=None, size=8.8):
    rh, hh = 16.0, 24.0
    total = sum(widths)
    x0 = (W - total) / 2.0
    if title:
        paragraph(title, size=size + 1.2, lead=1.6)
    need(hh + rh * (len(rows) + 1) + 8)
    ytop = Y[0]

    def put_cell(runs, cx, wd, ybase):
        wtot = runs_width(runs, size)
        draw_runs(runs, ybase, size, x_start=cx + wd / 2.0 + wtot / 2.0)

    cx = x0
    for hcell, wd in zip(reversed(headers), reversed(widths)):
        put_cell(cell_runs(hcell), cx, wd, ytop + 11)
        cx += wd
    page.draw_line((x0, ytop - 2), (x0 + total, ytop - 2), color=(0, 0, 0), width=0.9)
    page.draw_line((x0, ytop + hh - 8), (x0 + total, ytop + hh - 8), color=(0, 0, 0), width=0.9)
    yy = ytop + hh - 8
    for r in rows:
        yy += rh
        cx = x0
        for cell, wd in zip(reversed(r), reversed(widths)):
            put_cell(cell_runs(cell), cx, wd, yy + 11)
            cx += wd
        page.draw_line((x0, yy + rh - 4.5), (x0 + total, yy + rh - 4.5),
                       color=(0.6, 0.6, 0.6), width=0.4)
    cx = x0
    for wd in widths:
        page.draw_line((cx, ytop - 2), (cx, yy + rh - 4.5), color=(0.6, 0.6, 0.6), width=0.4)
        cx += wd
    page.draw_line((cx, ytop - 2), (cx, yy + rh - 4.5), color=(0.6, 0.6, 0.6), width=0.4)
    Y[0] = yy + rh + 8


def figure(path, caption, width=470.0):
    pix = pymupdf.Pixmap(path)
    hgt = width * pix.height / pix.width
    need(hgt + 60)
    page.insert_image(pymupdf.Rect((W - width) / 2, Y[0], (W + width) / 2, Y[0] + hgt),
                      filename=path)
    Y[0] += hgt + 14
    paragraph(caption, size=9.8, lead=1.6)
    Y[0] += 8


# ===================================================================== محتوا
Y[0] = 80.0
draw_runs([("fa", "اصلاحیه و الحاقیهٔ محاسبات میرایی فراصوت در ابررسانای هابارد با جفت‌شدگی اسپین‌مدار راشبا")],
          Y[0], 13.5)
Y[0] += 26
draw_runs([("fa", "نسخهٔ فارسی الحاقیهٔ مقالهٔ "), ("en", "arXiv:2607.02812")], Y[0], 11.5)
Y[0] += 26

heading("این سند چیست؟")
paragraph(("در مقالهٔ اصلی، نسبتِ توان جذب‌شده در حالت ابررسانا به توان جذب‌شده در حالت عادی، یعنی "
           "@@W/W_N@@، محاسبه و برای مجموعه‌های مختلفی از پارامترها رسم شده بود. آن محاسبه چند ایراد اساسی "
           "داشت. در فایل اصلاح‌شده ([[arXiv_2607.02812_corrected.pdf]]) متن اصلی مقاله دست‌نخورده باقی مانده "
           "و فقط هفت شکل نادرستِ صفحات ۱۱ تا ۱۵ با محاسبهٔ درست جایگزین شده‌اند و یک الحاقیهٔ شش‌صفحه‌ای "
           "(صفحات ۱۹ تا ۲۴) افزوده شده است. این سندِ فارسی، همان اصلاحیه‌ها را با زبانی ساده‌تر، با معادله‌های "
           "حروف‌چینی‌شده به سبک لاتک، و با همان جدول‌ها و شکل‌ها می‌آورد. در سراسر سند، @@\\hbar\\omega=0.7t@@ "
           "انرژی فوتون فراصوت و @@\\eta=+1@@ (مورد I، یعنی فراصوت) است."), size=10.5)
paragraph(("یک یادآوری فیزیکی کوتاه: فراصوت در یک ابررسانا با برانگیختن شبه‌ذره‌ها جذب می‌شود؛ یعنی یک "
           "شبه‌ذره از حالتی با انرژی @@E@@ به حالتی با انرژی @@E'=E+\\hbar\\omega@@ می‌رود (جذب فوتون). احتمال "
           "این گذار ضریبی به نام «ضریب هم‌دوسی» دارد که اثرِ ساختار ابررسانایی را وارد می‌کند. همهٔ بحث "
           "بر سر شکلِ درستِ همین ضریب و کران‌های درستِ انتگرال‌هاست."), size=10.5)

heading("به زبان ساده: هفت ایراد و هفت اصلاح")
paragraph(("بند ۱  ضریب هم‌دوسی.  ایراد: مقاله ضریب هم‌دوسی را با ریشه‌های نادرست و با جملهٔ نادرستِ شکاف "
           "ترکیب کرده بود و علامت @@\\xi@@ (فاصلهٔ انرژی از سطح شیمیایی) را نادیده می‌گرفت. اصلاح: شکل درستِ "
           "ضریب برای فراصوت، با حفظ علامت @@\\xi@@ و با پایستگی انرژی @@E'=E+\\hbar\\omega@@، معادلهٔ (۱) است؛ "
           "در این معادله @@\\bar{\\Delta}@@ شکاف وابسته به انرژی و جمع، روی هر دو شاخهٔ هلیسیتی راشبا انجام "
           "می‌شود:"), size=10.3)
eq_display(r"C_I(E,\xi;\,E',\xi') \;=\; \frac{1}{2}\left[\,1+\frac{\xi\,\xi'-\bar{\Delta}\,\bar{\Delta}'}"
           r"{E\,E'}\,\right]")
paragraph(("یک آزمون سادهٔ درستی: اگر به‌جای @@\\xi\\xi'@@ از @@|\\xi||\\xi'|@@ استفاده کنیم، در حدِ شکافِ صفر "
           "(@@\\Delta\\rightarrow 0@@) نسبت به بالای ۱ (تا ۱٫۹۷) می‌رسد؛ یعنی جذبِ ابررسانا از جذبِ حالتِ عادی "
           "بیشتر می‌شود بی‌آنکه کانال اضافه‌ای وجود داشته باشد. این غیرفیزیکی است و نشان می‌دهد شکل درست، "
           "همان معادلهٔ (۱) با علامتِ حفظ‌شده است."), size=10.3)
paragraph(("بند ۲  همهٔ عددها و منحنی‌ها.  ایراد: همهٔ مقدارهای @@W/W_N@@ که در مقاله آمده بود (منحنی‌های تخت "
           "بین ۰٫۹۰ و ۱٫۰۰) پیامد همان ایراد بند ۱ و ایراد بند ۵ است. اصلاح: منحنی‌های درست در همهٔ مجموعه‌ها "
           "کیفیتی مشابه دارند: در دمای پایین جذب نمایی کوچک است (شبه‌ذره‌ای نیست که برانگیخته شود)، با گرم‌شدن "
           "نسبت یکنوا بالا می‌آید، و در @@T\\rightarrow T_c^-@@ به مقدارِ زیرِ ۱ می‌رسد (جدول‌های ۱ و ۲ و شکل‌ها). "
           "برای مجموعهٔ پایه این حد برابر @@0.4024@@ و @@0.5086@@ و @@0.7616@@ به‌ترتیب برای @@V_{SO}=0@@ و "
           "@@0.4082@@ و @@0.8@@ است."), size=10.3)
paragraph(("بند ۳  حدِ نسبت در @@T_c@@.  ایراد: منحنی‌های مقاله دقیقاً در @@T_c@@ به ۱ می‌پریدند. اصلاح: در این "
           "مدل، ضریب هم‌دوسیِ هلیسیتی خاصیتِ خودِ نوارها (حالت عادی) است، نه خاصیتِ جفت‌شدگی؛ پس حتی در "
           "@@T_c@@ هم کانالِ پراکندگی کمتر از حالت عادی جذب می‌کند و نسبت در @@R(\\hbar\\omega,T_c)@@ (مقداری "
           "ساختار‌نواری، کوچک‌تر از ۱) اشباع می‌شود. کانالِ دوم، یعنی «شکستِ جفت»، در این کار تنها به‌عنوانِ "
           "چشم‌انداز بررسی می‌شود، نه نتیجهٔ کمّی. ضریب‌های هم‌دوسیِ دو کانال برای هر جفتِ حالتِ شبه‌ذره‌ای "
           "رابطهٔ (۲) را دقیقاً برآورده می‌کنند، پس وزنِ کلِ گذار بینِ دو کانال تقسیم می‌شود؛ اما چون دو کانال روی "
           "دو پوستهٔ متفاوت (@@E'-E=\\hbar\\omega@@ در برابر @@E+E'=\\hbar\\omega@@) و با وزن‌های گرماییِ متفاوت "
           "می‌نشینند، جمعِ آن‌ها با هیچ مرجعِ یگانهٔ حالتِ عادی برابر نمی‌شود. نتیجهٔ کمّیِ این کار فقط کانالِ "
           "پراکندگی است: کانالِ پراکندگی به‌تنهایی در حدِ @@T_c^-@@ مقدارهای @@0.400@@ و @@0.505@@ و @@0.761@@ را "
           "برمی‌گرداند؛ رسیدنِ @@W/W_N@@ به ۱ به کانالِ شکستِ جفت نیاز دارد و بیرون از دامنهٔ این نتیجهٔ کمّی است:"), size=10.3)
eq_display(r"C_{\mathrm{pair}}(E,E') \;=\; 1 - C_I(E,E')")
paragraph(("بند ۴  آستانهٔ شکستِ جفت.  ایراد: مقاله شرط @@2\\Delta_{eff}>\\hbar\\omega@@ را شرط رواییِ محاسبه "
           "گرفته بود. اصلاح: این شرط، آستانهٔ فیزیکیِ شکستِ جفت نیست؛ آستانهٔ واقعی معادلهٔ (۳) است، یعنی "
           "دو برابرِ کمینهٔ انرژی شبه‌ذره‌ای @@E_{min}@@. نکتهٔ مهم برای بحثِ روایی: درونِ دامنهٔ "
           "@@2\\Delta_{eff}>\\hbar\\omega@@ کانالِ شکستِ جفت عیناً صفر است؛ پس در همان دامنه‌ای که مقاله کار "
           "می‌کند، نتیجهٔ کانالِ پراکندگی به‌تنهایی کامل است و چیزی از قلم نیفتاده است:"), size=10.3)
eq_display(r"\hbar\omega > 2E_{min},\qquad E_{min}=\frac{|\Delta^0-b\,\mu|}{\sqrt{1+b^2}},\qquad "
           r"b=\frac{\Delta^s}{4t}")
paragraph(("بند ۵  مرجعِ حالتِ عادی و چگالی حالت‌ها.  ایراد: @@W@@ و @@W_N@@ باید با یک چگالی حالت‌ها "
           "@@N(\\varepsilon)@@ محاسبه شوند که نرمال‌سازیِ درست (@@\\int N\\,d\\varepsilon=2@@) و تکیه‌گاهِ دقیقِ "
           "نواری (@@|\\varepsilon|\\leq 4t\\sqrt{1+V_{SO}^2/2t^2}@@) داشته باشد؛ کد اصلی نرمال‌سازیِ نادرست "
           "(@@0.643@@) و نشتِ گاوسی بالای لبهٔ نوار داشت و مقدارهای منفی را با قدرمطلق پنهان می‌کرد. اصلاح: "
           "هر دو انتگرال با یک چگالیِ حالتِ درست و یکسان محاسبه شدند."), size=10.3)
paragraph(("بند ۶  دو نادرستیِ پارامتری.  عنوان گروه B نوشته بود @@\\Delta_t=4@@، ولی عددهای چاپ‌شدهٔ خودِ "
           "مقاله تنها با @@\\Delta_t=4.5@@ بازتولید می‌شوند (غلط تایپی). پنل‌های کمکیِ گروه D هم عددهایی را "
           "نشان می‌دادند که متعلق به گروه D نبودند (عیناً سطر @@U=-4t@@ گروه A بودند). هر دو بازتولید و "
           "اصلاح شدند."), size=10.3)
paragraph(("بند ۷  مجموعهٔ @@U=-75t@@.  برای این مجموعه، حلِ بیشینه‌شکافِ معادله‌های مرحلهٔ ۱ (همان قاعدهٔ "
           "انتخابِ حل‌گرِ خودِ مقاله) یک ریشهٔ جفت‌شدگیِ قوی می‌دهد با سطح شیمیایی @@\\mu\\approx 35t@@، یعنی "
           "بسیار بیرون از پهنای نوار (@@T_c\\approx 12.7t@@). چنین حلّی در این مدلِ تنگ‌بستن فیزیکی نیست؛ "
           "منحنیِ اصلاح‌شده‌اش همان‌گونه که هست رسم شده، ولی نباید بر آن تکیهٔ فیزیکی کرد."), size=10.3)

heading("صورت‌بندی و معادله‌های اصلاح‌شده")
paragraph(("مدل: شبکهٔ مربعی با پرش @@t@@، جفت‌شدگی اسپین‌مدار راشبا @@V_{SO}@@ و برهم‌کنش هابارد @@U@@. "
           "انرژی نوارها معادلهٔ (۴) و شکافِ وابسته به انرژی (فرم HHM) و انرژی شبه‌ذره‌ای معادلهٔ (۵) است:"),
          size=10.3)
eq_display(r"\varepsilon_{\pm}(k)=-2t(\cos k_x+\cos k_y)\mp 2V_{SO}|\sin k|")
eq_display(r"\bar{\Delta}(\varepsilon)=\Delta^0-\frac{\Delta^s}{4t}\,\varepsilon,\qquad "
           r"E_{\sigma}(\varepsilon)=\sqrt{(\varepsilon-\mu)^2+\bar{\Delta}(\varepsilon)^2}")
paragraph(("با معرفیِ پارامترهای معادلهٔ (۶)، رابطهٔ انرژی نوار با انرژی شبه‌ذره‌ای برای هر شاخه "
           "(@@\\sigma=\\pm 1@@) دقیقاً وارون می‌شود و measure انتگرال ساده می‌شود؛ این کار تکینگیِ "
           "@@1/\\sqrt{E^2-E_{min}^2}@@ را دقیقاً حذف می‌کند و انتگرال‌ها را پایدار می‌سازد:"), size=10.3)
eq_display(r"b=\frac{\Delta^s}{4t},\quad \kappa=\sqrt{1+b^2},\quad "
           r"\varepsilon^{*}=\frac{\mu+b\Delta^0}{1+b^2},\quad \ell=\sqrt{E^2-E_{min}^2}")
eq_display(r"\varepsilon_{\sigma}(E)=\varepsilon^{*}+\sigma\frac{\ell}{\kappa},\qquad "
           r"\nu_{\sigma}(E)\,dE=N(\varepsilon_{\sigma}(E))\,\frac{d\ell}{\kappa}", numbered=False)
paragraph(("کانالِ پراکندگی (گذارِ یک شبه‌ذره از @@E@@ به @@E+\\hbar\\omega@@) با ضریب هم‌دوسیِ معادلهٔ (۱) "
           "می‌شود معادلهٔ (۷)؛ کرانِ بالای هر جفتِ شاخه از تکیه‌گاهِ نوار می‌آید: "
           "@@E\\leq\\min(E_{max}(\\sigma),\\,E_{max}(\\sigma')-\\hbar\\omega)@@. مرجعِ حالتِ عادی معادلهٔ (۸) "
           "و کانالِ شکستِ جفت معادلهٔ (۹) است:"), size=10.3)
eq_display(r"W=\sum_{\sigma,\sigma'}\int dE\ \nu_{\sigma}(E)\,\nu_{\sigma'}(E+\hbar\omega)\ "
           r"C_I\ \left[f(E)-f(E+\hbar\omega)\right]")
eq_display(r"W_N=\int dE\ N(E)\,N(E+\hbar\omega)\left[f_N(E)-f_N(E+\hbar\omega)\right]")
eq_display(r"W_p=\sum_{\sigma,\sigma'}\int dE\ \nu_{\sigma}(E)\,\nu_{\sigma'}(\hbar\omega-E)\ "
           r"(1-C_I)\ \left[1-f(E)-f(\hbar\omega-E)\right]")
paragraph(("در @@\\hbar\\omega=0.7t@@ دو جفتِ شاخه @@(+1,+1)@@ و @@(-1,+1)@@ عیناً صفرند، زیرا بیشینهٔ انرژی "
           "شاخهٔ بالایی از انرژی فوتون کمتر است. حدِ @@\\Delta\\rightarrow 0@@ نسبتِ پراکندگی هم فرم بستهٔ "
           "معادلهٔ (۱۰) را دارد (با دقت @@3\\times10^{-6}@@ تأیید شده است):"), size=10.3)
eq_display(r"R(\hbar\omega,T)=\frac{2}{\hbar\omega}\int_0^{\infty}\left[f(E)-f(E+\hbar\omega)"
           r"\right]dE")

heading("جدول‌های نتیجهٔ اصلاح‌شده")
paragraph(("جدول ۱ همهٔ مجموعه‌های گروه‌های A تا D را خلاصه می‌کند. ستونِ «کمینهٔ نسبت» نشان می‌دهد در دمای "
           "پایین جذب عملاً صفر است؛ ستونِ آخر، مقدارِ اشباع در @@T\\rightarrow T_c^-@@ را نشان می‌دهد که در همهٔ "
           "سطرها زیر ۱ است:"), size=10.3)
table([("fa", "مجموعه"), ("eq", "T_c/t"), ("eq", "\\Delta_{eff}(0)/t"), ("fa", "کمینهٔ نسبت"),
       ("eq", "T\\rightarrow T_c^-")],
      [[("en", s), ("en", a), ("en", b), ("en", "0.000"), ("en", c)]
       for s, a, b, c in [
           ("A:  U=-t", "0.0123", "0.0214", "0.0198"), ("A:  U=-2t", "0.1818", "0.3206", "0.2949"),
           ("A:  U=-3t", "0.4388", "0.7969", "0.5567"), ("A:  U=-4t", "0.7184", "1.3273", "0.6731"),
           ("B:  U=75t", "0.5009", "0.6272", "0.6717"), ("B:  U=115t", "0.1399", "0.2840", "0.1837"),
           ("B:  U=125t", "0.0861", "0.2070", "0.1063"), ("C:  Δ_t=3.5", "0.4094", "0.6057", "0.4652"),
           ("C:  Δ_t=4.0", "0.6720", "0.8414", "0.7930"), ("C:  Δ_t=4.5", "0.9536", "1.0934", "0.9890"),
           ("D:  V_SO=0", "0.2861", "0.4444", "0.3999"), ("D:  V_SO=0.4", "0.3770", "0.5447", "0.5019"),
           ("D:  V_SO=0.8", "0.6494", "0.8542", "0.7602")]],
      widths=[96, 62, 84, 84, 84],
      title="جدول ۱: گروه‌های A تا D، کانال پراکندگی با ضریب درستِ معادلهٔ (۱)، @@\\hbar\\omega=0.7t@@.")
table([("eq", "V_{SO}"), ("eq", "T_c/t"), ("eq", "\\Delta_{eff}(0)/t"), [("fa", "حدِ نسبت در "), ("eq", "T_c")],
       ("fa", "مقدار ساختار‌نواری")],
      [[("en", "0.000"), ("en", "0.2861"), ("en", "0.4444"), ("en", "0.40241"), ("en", "0.49867")],
       [("en", "0.4082"), ("en", "0.3808"), ("en", "0.5489"), ("en", "0.50856"), ("en", "0.59355")],
       [("en", "0.800"), ("en", "0.6494"), ("en", "0.8542"), ("en", "0.76163"), ("en", "0.74264")]],
      widths=[66, 66, 88, 104, 104],
      title="جدول ۲: مجموعهٔ پایه (@@U=95t@@، @@\\Delta_t=4.5t@@، @@n=1.885@@)؛ ستون آخر مقدارِ تحلیلیِ "
            "@@R(0.7,T_c)@@ با چگالیِ حالتِ تخت است که حدِ عددی را تأیید می‌کند.")
table([("fa", "مجموعه"), [("fa", "گسترهٔ نسبت از "), ("eq", "V_{SO}=0"), ("fa", " تا "), ("eq", "0.8")]],
      [[("en", "U=95t"), ("en", "0.253   →   0.435")],
       [("en", "U=125t"), ("en", "0.069   →   0.253")],
       [("en", "U=-75t"), [("en", "0.489   →   0.496"), ("fa", " (بند ۷)")]]],
      widths=[100, 220],
      title="جدول ۳: جاروب‌های @@V_{SO}@@ در @@T/T_c=0.9@@؛ نسبت با @@V_{SO}@@ و با دما بالا می‌آید ولی هرگز "
            "به ۱ نمی‌رسد.")
paragraph(("جمع‌بندی جدول‌ها: (i) در هر سطر و هر دما @@W/W_N<1@@؛ (ii) حدِ @@T_c@@ برابر مقدارِ ساختار‌نواریِ "
           "@@R(\\hbar\\omega,T_c)@@ است نه ۱ (بند ۳)؛ (iii) درونِ دامنهٔ رواییِ مقاله کانالِ شکستِ جفت صفر است "
           "(بند ۴)؛ برای مجموعهٔ پایه آستانهٔ واقعیِ (۳) در @@T/T_c=0.820@@ و @@0.912@@ و @@0.975@@ است، "
           "حال‌آنکه برشِ مقاله (@@2\\Delta_{eff}=\\hbar\\omega@@) در ۰٫۷۴۰ و ۰٫۸۴۰ و ۰٫۹۳۶ می‌افتد؛ (iv) بالای "
           "آستانه، کانالِ شکستِ جفت (۹) تنها چشم‌انداز است و در نتیجهٔ کمّی وارد نمی‌شود (بند ۳)."), size=10.3)

heading("اعتبارسنجی: از کجا می‌دانیم محاسبه درست است؟")
paragraph(("۱) انتگرال‌گیری با دو روش مستقل (کوادراتور تطبیقی و قانون نقطهٔ میانی) تا اختلاف نسبیِ "
           "@@2\\times10^{-4}@@ توافق دارد. ۲) در حدِ شکافِ صفر، نتیجه با فرم بستهٔ (۱۰) تا "
           "@@3\\times10^{-6}@@ برابر است و رابطهٔ (۲) برای هر جفتِ حالت دقیقاً برقرار است (در حدِ شکافِ صفر، جمعِ دو "
           "کانال برابر مرجعِ حالت عادی می‌شود — آزمون T1). ۳) یک شمارشِ مستقلِ "
           "شبکه‌ای روی جفت‌های @@(k,k')@@ (بی‌آنکه از هیچ‌کدام از انتگرال‌های بالا استفاده کند) با نتیجهٔ "
           "روش اصلی تا ۰٫۷ تا ۱٫۳ درصد در نزدیکی @@T_c@@ توافق دارد. ۴) ضریب هم‌دوسیِ (۱) در سراسر دامنهٔ "
           "انتگرال بین ۰ و ۱ است، همان‌گونه که باید."), size=10.3)

heading("بازتولید محاسبه‌ها")
paragraph(("مرحلهٔ ۱ (یافتنِ @@T_c@@ و شکاف و سطح شیمیایی هر مجموعه): اسکریپت [[paper_fix/stage1.py]] که "
           "رونوشتِ لفظ‌به‌لفظِ حل‌گرِ مخزنِ خودِ مقاله است و همهٔ عددهای چاپ‌شدهٔ مقاله را تا چهار رقم اعشار "
           "بازتولید می‌کند. مرحلهٔ ۲ (نسبتِ @@W/W_N@@): خط لولهٔ اعتبارسنجی‌شدهٔ [[wwn_eta1/model.py]] و "
           "[[wwn_eta1/energy_space.py]] با آزمون‌های بالا در [[wwn_eta1/validate.py]]؛ اجراهای گروه‌ها و "
           "جاروب‌ها: [[paper_fix/run_groups.py]] و [[paper_fix/run_scans.py]]؛ شکل‌ها: "
           "[[paper_fix/make_paper_figs.py]]؛ فایل اصلاح‌شدهٔ انگلیسی: [[paper_fix/build_pdf.py]] و همین سندِ "
           "فارسی: [[paper_fix/build_pdf_fa.py]]. گزارش کاملِ استخراج و اعتبارسنجی هم در "
           "[[wwn_eta1/REPORT.md]] آمده است."), size=10.3)

heading("شکل‌های مجموعهٔ پایه")
figure(ETA1 + "fig_WWN_vs_T.png",
       ("شکل ۱: نسبتِ اصلاح‌شدهٔ @@W/W_N@@ بر حسب @@T/T_c@@ برای سه مقدارِ @@V_{SO}@@ در مجموعهٔ پایه. "
        "منحنی‌های رنگی کم‌رنگ‌تر، شکل‌های نادرستِ ضریب هم‌دوسی‌اند (گزینه‌های (b) و (c) و (d)) که یا از ۱ "
        "عبور می‌کنند یا حدِ نادرست دارند؛ تنها منحنی‌های پررنگ (گزینهٔ (a)، یعنی معادلهٔ (۱)) هم زیر ۱ "
        "می‌مانند و هم حدِ درست دارند. کانالِ شکستِ جفت رسم نشده است (چشم‌اندازِ بدونِ ادّعای کمّی، بند ۳)."))
figure(ETA1 + "fig_branch_decomposition.png",
       ("شکل ۲: سهمِ هر جفتِ شاخهٔ هلیسیتی @@(\\sigma,\\sigma')@@ از جمعِ پراکندگی. کانالِ درون‌نواریِ "
        "@@(-1,-1)@@ غالب است و دو کانالِ دیگر در @@\\hbar\\omega=0.7t@@ دقیقاً صفرند؛ این رفتارِ صفر‌شدن، "
        "یکی از بررسی‌های درستیِ کران‌های انتگرال است."))
figure(ETA1 + "fig_limits.png",
       ("شکل ۳: دو حدِ تحلیلی. راست: در حدِ شکافِ صفر، نتیجهٔ عددی روی منحنیِ تحلیلیِ (۱۰) می‌نشیند و با "
        "گفتهٔ کتاب‌درسیِ @@2f(\\Delta)@@ مقایسه شده است (تفاوت، به تفاوتِ نقطهٔ ارزیابیِ ضریب برمی‌گردد، نه "
        "به ایرادِ کد). چپ: اشباعِ نسبت در @@R(\\hbar\\omega,T_c)@@ و بازگشتِ دقیق به ۱ پس از افزودنِ کانالِ "
        "شکستِ جفت بالای آستانهٔ (۳)."))

for i, pg in enumerate(doc):
    n = fad(i + 1)
    w = w_fa(n, 9)
    pg.insert_text(((W - w) / 2, H - 28), fa_vis(n), fontname="fa", fontfile=FA_TTF,
                   fontsize=9, color=(0.3, 0.3, 0.3))

doc.set_metadata({"title": "الحاقیهٔ اصلاحی arXiv:2607.02812 (نسخهٔ فارسی)", "author": "paper_fix"})
doc.save(os.path.join(HERE, "results", "corrections_fa.pdf"), garbage=4, deflate=True)
print("saved corrections_fa.pdf, pages:", doc.page_count, "equations:", EQNO[0])
