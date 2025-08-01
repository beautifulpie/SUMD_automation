#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CSV 파일 검증 및 전처리 도구

CSV 파일의 형식과 내용을 검증하고 필요시 정제합니다.
"""

import pandas as pd
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

class CSVValidator:
    """CSV 파일 검증 클래스"""
    
    def __init__(self):
        self.required_columns = ['PDB', 'Peptide Chain', 'Receptor Chain']
        self.valid_pdb_pattern = re.compile(r'^[0-9][a-zA-Z0-9]{3}$')
        self.valid_chain_pattern = re.compile(r'^[A-Za-z0-9]$')
    
    def validate_csv(self, csv_file: str) -> Dict:
        """CSV 파일 전체 검증"""
        
        result = {
            'valid': False,
            'errors': [],
            'warnings': [],
            'stats': {},
            'cleaned_entries': 0
        }
        
        try:
            # CSV 로드
            df = pd.read_csv(csv_file)
            result['stats']['total_rows'] = len(df)
            
            print(f"CSV 파일 로드: {len(df)} 행")
            
            # 필수 컬럼 확인
            missing_columns = [col for col in self.required_columns if col not in df.columns]
            if missing_columns:
                result['errors'].append(f"필수 컬럼 누락: {missing_columns}")
                return result
            
            # 각 행 검증
            valid_rows = []
            invalid_count = 0
            
            for idx, row in df.iterrows():
                validation = self.validate_row(row, idx)
                
                if validation['valid']:
                    valid_rows.append(idx)
                else:
                    invalid_count += 1
                    for error in validation['errors']:
                        result['errors'].append(f"행 {idx+1}: {error}")
                
                for warning in validation['warnings']:
                    result['warnings'].append(f"행 {idx+1}: {warning}")
            
            result['stats']['valid_rows'] = len(valid_rows)
            result['stats']['invalid_rows'] = invalid_count
            result['stats']['validity_rate'] = len(valid_rows) / len(df) * 100
            
            # 전체 유효성
            if len(valid_rows) > 0:
                result['valid'] = True
                result['cleaned_entries'] = len(valid_rows)
            
            # 추가 통계
            if result['valid']:
                result['stats'].update(self._compute_additional_stats(df.iloc[valid_rows]))
            
            return result
            
        except Exception as e:
            result['errors'].append(f"CSV 파일 로드 실패: {e}")
            return result
    
    def validate_row(self, row: pd.Series, row_idx: int) -> Dict:
        """개별 행 검증"""
        
        validation = {
            'valid': True,
            'errors': [],
            'warnings': []
        }
        
        # PDB ID 검증
        pdb_id = str(row.get('PDB', '')).strip().lower()
        if not pdb_id:
            validation['errors'].append("PDB ID 누락")
            validation['valid'] = False
        elif not self.valid_pdb_pattern.match(pdb_id):
            validation['errors'].append(f"유효하지 않은 PDB ID: {pdb_id}")
            validation['valid'] = False
        
        # Peptide Chain 검증
        peptide_chain = str(row.get('Peptide Chain', '')).strip()
        if not peptide_chain or peptide_chain.lower() in ['nan', 'none']:
            validation['errors'].append("Peptide Chain 누락")
            validation['valid'] = False
        elif not self.valid_chain_pattern.match(peptide_chain):
            validation['errors'].append(f"유효하지 않은 Peptide Chain: {peptide_chain}")
            validation['valid'] = False
        
        # Receptor Chain 검증
        receptor_chain = str(row.get('Receptor Chain', '')).strip()
        if not receptor_chain or receptor_chain.lower() in ['nan', 'none']:
            validation['errors'].append("Receptor Chain 누락")
            validation['valid'] = False
        elif not self.valid_chain_pattern.match(receptor_chain):
            validation['errors'].append(f"유효하지 않은 Receptor Chain: {receptor_chain}")
            validation['valid'] = False
        
        # Chain 중복 확인
        if validation['valid'] and peptide_chain == receptor_chain:
            validation['warnings'].append("Peptide와 Receptor Chain이 동일함")
        
        return validation
    
    def clean_csv(self, input_csv: str, output_csv: str) -> Dict:
        """CSV 파일 정제"""
        
        print(f"CSV 정제 시작: {input_csv} -> {output_csv}")
        
        try:
            df = pd.read_csv(input_csv)
            original_count = len(df)
            
            # 유효한 행들만 필터링
            valid_indices = []
            
            for idx, row in df.iterrows():
                validation = self.validate_row(row, idx)
                if validation['valid']:
                    valid_indices.append(idx)
            
            # 정제된 DataFrame 생성
            cleaned_df = df.iloc[valid_indices].copy()
            
            # PDB ID 정규화 (소문자)
            cleaned_df['PDB'] = cleaned_df['PDB'].str.lower().str.strip()
            
            # Chain ID 정규화 (대문자)
            cleaned_df['Peptide Chain'] = cleaned_df['Peptide Chain'].str.upper().str.strip()
            cleaned_df['Receptor Chain'] = cleaned_df['Receptor Chain'].str.upper().str.strip()
            
            # 중복 제거 (PDB + Chain 조합 기준)
            cleaned_df = cleaned_df.drop_duplicates(
                subset=['PDB', 'Peptide Chain', 'Receptor Chain'],
                keep='first'
            )
            
            # 정제된 파일 저장
            cleaned_df.to_csv(output_csv, index=False)
            
            result = {
                'success': True,
                'original_count': original_count,
                'cleaned_count': len(cleaned_df),
                'removed_count': original_count - len(cleaned_df),
                'removal_rate': (original_count - len(cleaned_df)) / original_count * 100
            }
            
            print(f"정제 완료: {original_count} -> {len(cleaned_df)} 항목")
            print(f"제거율: {result['removal_rate']:.1f}%")
            
            return result
            
        except Exception as e:
            return {
                'success': False,
                'error': str(e)
            }
    
    def _compute_additional_stats(self, df: pd.DataFrame) -> Dict:
        """추가 통계 계산"""
        
        stats = {}
        
        # 고유 PDB 수
        stats['unique_pdbs'] = int(df['PDB'].nunique())
        
        # Chain 통계
        stats['unique_peptide_chains'] = int(df['Peptide Chain'].nunique())
        stats['unique_receptor_chains'] = int(df['Receptor Chain'].nunique())
        
        # Chain 분포 (JSON 직렬화 가능하도록 변환)
        peptide_chain_counts = df['Peptide Chain'].value_counts().to_dict()
        receptor_chain_counts = df['Receptor Chain'].value_counts().to_dict()
        
        # 모든 값을 int로 변환 (numpy int64 등을 방지)
        stats['peptide_chain_distribution'] = {str(k): int(v) for k, v in peptide_chain_counts.items()}
        stats['receptor_chain_distribution'] = {str(k): int(v) for k, v in receptor_chain_counts.items()}
        
        # 가장 일반적인 조합 (tuple 키를 string으로 변환)
        try:
            common_combinations = df.groupby(['Peptide Chain', 'Receptor Chain']).size().sort_values(ascending=False)
            combinations_dict = {}
            
            for (peptide, receptor), count in common_combinations.head(10).items():
                key = f"{peptide}-{receptor}"  # tuple을 string으로 변환
                combinations_dict[key] = int(count)
            
            stats['common_chain_combinations'] = combinations_dict
        except Exception as e:
            # 조합 계산 실패시 빈 딕셔너리
            stats['common_chain_combinations'] = {}
        
        return stats


def main():
    """메인 함수"""
    import argparse
    
    parser = argparse.ArgumentParser(description="CSV 파일 검증 및 정제")
    parser.add_argument("--input", required=True, help="입력 CSV 파일")
    parser.add_argument("--validate", action="store_true", help="검증만 수행")
    parser.add_argument("--clean", help="정제된 CSV 출력 파일")
    parser.add_argument("--report", help="검증 보고서 출력 파일")
    
    args = parser.parse_args()
    
    if not Path(args.input).exists():
        print(f"❌ 입력 파일을 찾을 수 없습니다: {args.input}")
        return 1
    
    validator = CSVValidator()
    
    # 검증 수행
    print("=== CSV 검증 시작 ===")
    validation_result = validator.validate_csv(args.input)
    
    # 검증 결과 출력
    print(f"\n=== 검증 결과 ===")
    print(f"전체 유효성: {'✅ 유효' if validation_result['valid'] else '❌ 무효'}")
    print(f"총 행 수: {validation_result['stats'].get('total_rows', 0)}")
    print(f"유효한 행: {validation_result['stats'].get('valid_rows', 0)}")
    print(f"무효한 행: {validation_result['stats'].get('invalid_rows', 0)}")
    print(f"유효율: {validation_result['stats'].get('validity_rate', 0):.1f}%")
    
    if validation_result['stats'].get('unique_pdbs'):
        print(f"고유 PDB 수: {validation_result['stats']['unique_pdbs']}")
    
    # 오류 출력
    if validation_result['errors']:
        print(f"\n❌ 오류 ({len(validation_result['errors'])}개):")
        for error in validation_result['errors'][:10]:  # 처음 10개만
            print(f"  - {error}")
        if len(validation_result['errors']) > 10:
            print(f"  ... 및 {len(validation_result['errors']) - 10}개 더")
    
    # 경고 출력
    if validation_result['warnings']:
        print(f"\n⚠️  경고 ({len(validation_result['warnings'])}개):")
        for warning in validation_result['warnings'][:5]:  # 처음 5개만
            print(f"  - {warning}")
        if len(validation_result['warnings']) > 5:
            print(f"  ... 및 {len(validation_result['warnings']) - 5}개 더")
    
    # 정제 수행
    if args.clean and validation_result['valid']:
        print(f"\n=== CSV 정제 시작 ===")
        clean_result = validator.clean_csv(args.input, args.clean)
        
        if clean_result['success']:
            print(f"✅ 정제 완료: {args.clean}")
            print(f"원본: {clean_result['original_count']} 항목")
            print(f"정제 후: {clean_result['cleaned_count']} 항목")
            print(f"제거: {clean_result['removed_count']} 항목 ({clean_result['removal_rate']:.1f}%)")
        else:
            print(f"❌ 정제 실패: {clean_result['error']}")
    
    # 보고서 저장
    if args.report:
        import json
        with open(args.report, 'w', encoding='utf-8') as f:
            json.dump(validation_result, f, indent=2, ensure_ascii=False, default=str)
        print(f"📋 검증 보고서 저장: {args.report}")
    
    return 0 if validation_result['valid'] else 1


if __name__ == "__main__":
    exit(main())