#!/bin/sh
shopt -s nocasematch

sep=""
format="standard"
rec=1
help=0
hours=0
minutes=0
seconds=0
years=0

while getopts "Hbehmsf:t:y" opt; do
    case $opt in
        H)
            help=1
            ;;
        b)
            rec=2
            ;;
        e)
            rec=1
            ;;
        h)
            hours=1
            ;;
        m)
            hours=1
            minutes=1
            ;;
        s)
            hours=1
            minutes=1
            seconds=1
            ;;
        f)
            format="$OPTARG"
            ;;
        t)
            sep="$OPTARG"
            ;;
        y)
            years=1
            ;;
    esac
done

# --- help output ---
if (( help )); then
    cat <<END
time_stamp.sh [ -Hbehms -f format -t separator ]

  -H help (no execution)
  -b beginning date
  -e ending date (default)
  -h hours
  -m hours,minutes
  -s hours,minutes,seconds
  -f format=standard(default),european,digital,days
  -t separator (default=blank)
  -y year

END
    exit
fi

# --- check format ---
if [[ $format != "standard" && $format != "european" && $format != "digital" && $format != "days" ]]; then
    echo "ERROR: Invalid format"
    exit 4
fi

hsep=$sep
if [[ $format == "standard" || $format == "european" || $format == "days" ]]; then
    hsep="h"
fi

if [[ -e time_stamp.out ]]; then
    time_stamp=($(tail -n "$rec" time_stamp.out))

    month_name=$(echo "$time_stamp" | awk '{print tolower($7)}')
    month_num=$(printf "%.2d" "${time_stamp[1]}")

    # ---- day can have more than 2 digits ----
    day_num=$(printf "%.2d" "${time_stamp[2]}")
    if [[ $month_name == "day" && $format == "standard" ]]; then
        day_num=$(printf "%.4d" "${time_stamp[2]}")
    fi

    # ---- hours,min,sec can have only 2 digits ----
    hour_num=$(printf "%.2d" "${time_stamp[3]}")
    min_num=$(printf "%.2d" "${time_stamp[4]}")
    sec_num=$(printf "%.2d" "${time_stamp[5]}")

    # ---- pad ISO years to 4 digits ----
    year=""
    if [[ $format == "digital" ]]; then
        year=$(printf "%.4d" "${time_stamp[0]}") # will work even if year > 9999
    fi

    # ---- pad days to 5 digits for no_calendar models ----
    if [[ $format == "days" ]]; then
        day_num=$(printf "%.5d" "${time_stamp[2]}")
    fi

    # ---- create date label ----
    date_name=""

    if [[ $format == "standard" ]]; then
        if [[ $month_name != "day" ]]; then
            date_name="${time_stamp[0]}"
        fi
        date_name="$date_name$month_name$day_num"
    elif [[ $format == "european" ]]; then
        date_name="$day_num$month_name${time_stamp[0]}"
    elif [[ $format == "digital" ]]; then
        date_name="$year$sep$month_num$sep$day_num"
    elif [[ $format == "days" ]]; then
        date_name="$day_num"
    fi

    if (( hours )); then
        date_name="$date_name$hsep$hour_num"
    fi
    if (( minutes )); then
        date_name="$date_name$sep$min_num"
    fi
    if (( seconds )); then
        date_name="$date_name$sep$sec_num"
    fi
    if (( years )); then
        date_name="$year"
    fi
else
    # --- dummy values ---
    month_name="xxx"
    date_name="no_time_stamp"
fi

echo "$date_name"