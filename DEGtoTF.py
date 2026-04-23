#!/usr/bin/env python3

import argparse
import pandas as pd
import os
import sys


def read_table(path):
    """Чтение таблицы (xlsx или csv)"""
    if not os.path.exists(path):
        sys.exit(f"Ошибка: файл не найден -> {path}")

    if path.endswith(".xlsx") or path.endswith(".xls"):
        return pd.read_excel(path, engine="openpyxl")
    elif path.endswith(".csv"):
        return pd.read_csv(path)
    else:
        sys.exit(f"Ошибка: неподдерживаемый формат файла -> {path}")


def main():
    parser = argparse.ArgumentParser(description="Mapping DEG to TF")

    parser.add_argument("-TF", required=True, help="TF table")
    parser.add_argument("-DEG", action="append", required=True, help="DEG tables (multiple allowed)")
    parser.add_argument("-OUT", required=True, help="Output file (.xlsx)")
    parser.add_argument("-LogFC", type=float, required=True, help="LogFC threshold (absolute)")
    parser.add_argument("-FDR", type=float, required=True, help="FDR threshold")

    args = parser.parse_args()

    # === 1. Чтение TF ===
    tf_df = read_table(args.TF)

    if "Locus Tag" not in tf_df.columns:
        sys.exit("Ошибка: в TF таблице нет столбца 'Locus Tag'")

    # === 2. Подготовка результирующего DataFrame ===
    result_df = tf_df.copy()

    # Счётчик совпадений
    matched_tf = set()

    # === 3. Обработка DEG таблиц ===
    for i, deg_path in enumerate(args.DEG, start=1):
        deg_df = read_table(deg_path)

        # Проверка столбцов
        for col in ["locus", "logFC", "FDR"]:
            if col not in deg_df.columns:
                sys.exit(f"Ошибка: в {deg_path} нет столбца '{col}'")

        # Оставляем только нужные столбцы
        deg_sub = deg_df[["locus", "logFC", "FDR"]].copy()

        # Переименовываем
        deg_sub = deg_sub.rename(columns={
            "locus": "Locus Tag",
            "logFC": f"LogFC{i}",
            "FDR": f"FDR{i}"
        })

        # Merge
        result_df = result_df.merge(deg_sub, on="Locus Tag", how="left")

        # Обновляем список совпадений
        matched = result_df[~result_df[f"LogFC{i}"].isna()]["Locus Tag"]
        matched_tf.update(matched)

    print(f"Сопоставление завершено: найдено совпадений для {len(matched_tf)} TF")

    # === 4. Фильтрация ===
    filtered_df = result_df.copy()

    for i in range(1, len(args.DEG) + 1):
        logfc_col = f"LogFC{i}"
        fdr_col = f"FDR{i}"

        # Условия:
        # 1. значение существует (не NaN)
        # 2. FDR <= threshold
        # 3. abs(LogFC) >= threshold
        condition = (
            filtered_df[logfc_col].notna() &
            filtered_df[fdr_col].notna() &
            (filtered_df[fdr_col] <= args.FDR) &
            (filtered_df[logfc_col].abs() >= args.LogFC)
        )

        filtered_df = filtered_df[condition]

    print(f"После фильтрации осталось {len(filtered_df)} TF со статистически значимыми значениями")

    # === 5. Сохранение ===
    try:
        filtered_df.to_excel(args.OUT, index=False, engine="openpyxl")
    except Exception as e:
        sys.exit(f"Ошибка при сохранении файла: {e}")

    print(f"Результат сохранён в файл: {args.OUT}")


if __name__ == "__main__":
    main()
