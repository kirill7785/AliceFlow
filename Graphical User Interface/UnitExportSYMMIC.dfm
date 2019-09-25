object Formexport2SYMMIC: TFormexport2SYMMIC
  Left = 322
  Top = 124
  Width = 476
  Height = 147
  Caption = 'export to SYMMIC'
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object lbl1: TLabel
    Left = 16
    Top = 8
    Width = 254
    Height = 13
    Caption = 'Enter the file name without an extension for SYMMIC'
  end
  object lbl2: TLabel
    Left = 16
    Top = 40
    Width = 143
    Height = 13
    Caption = 'Enter the name of the header'
  end
  object lblnx: TLabel
    Left = 16
    Top = 80
    Width = 12
    Height = 13
    Caption = 'nx'
  end
  object lblny: TLabel
    Left = 120
    Top = 80
    Width = 12
    Height = 13
    Caption = 'ny'
  end
  object lblnz: TLabel
    Left = 224
    Top = 80
    Width = 11
    Height = 13
    Caption = 'nz'
  end
  object lbltol: TLabel
    Left = 304
    Top = 48
    Width = 45
    Height = 13
    Caption = 'tolerance'
  end
  object edtname: TEdit
    Left = 280
    Top = 8
    Width = 169
    Height = 21
    TabOrder = 0
  end
  object edttitle: TEdit
    Left = 168
    Top = 40
    Width = 121
    Height = 21
    TabOrder = 1
  end
  object cbbnx: TComboBox
    Left = 40
    Top = 72
    Width = 57
    Height = 21
    ItemHeight = 13
    ItemIndex = 0
    TabOrder = 2
    Text = '5'
    Items.Strings = (
      '5'
      '10'
      '15'
      '20'
      '25')
  end
  object cbbny: TComboBox
    Left = 144
    Top = 72
    Width = 57
    Height = 21
    ItemHeight = 13
    ItemIndex = 0
    TabOrder = 3
    Text = '5'
    Items.Strings = (
      '5'
      '10'
      '15'
      '20'
      '25')
  end
  object cbbnz: TComboBox
    Left = 248
    Top = 72
    Width = 57
    Height = 21
    ItemHeight = 13
    ItemIndex = 0
    TabOrder = 4
    Text = '5'
    Items.Strings = (
      '5'
      '10'
      '15'
      '20'
      '25')
  end
  object cbbtol: TComboBox
    Left = 360
    Top = 40
    Width = 89
    Height = 21
    ItemHeight = 13
    ItemIndex = 0
    TabOrder = 5
    Text = '5%'
    Items.Strings = (
      '5%'
      '10%'
      '15%'
      '20%'
      '25%'
      '30%')
  end
  object btnApply: TButton
    Left = 344
    Top = 72
    Width = 75
    Height = 25
    Caption = 'Apply'
    TabOrder = 6
    OnClick = btnApplyClick
  end
end
