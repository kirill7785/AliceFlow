object FormViewFactors: TFormViewFactors
  Left = 0
  Top = 0
  AutoSize = True
  Caption = 'ViewFactors'
  ClientHeight = 497
  ClientWidth = 473
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object GroupBox1: TGroupBox
    Left = 0
    Top = 176
    Width = 473
    Height = 113
    Caption = 'View Factors'
    TabOrder = 0
    object Label1: TLabel
      Left = 16
      Top = 72
      Width = 31
      Height = 13
      Caption = 'Label1'
    end
    object Label2: TLabel
      Left = 88
      Top = 72
      Width = 31
      Height = 13
      Caption = 'Label2'
    end
    object Label3: TLabel
      Left = 168
      Top = 72
      Width = 31
      Height = 13
      Caption = 'Label3'
    end
    object Label4: TLabel
      Left = 240
      Top = 72
      Width = 31
      Height = 13
      Caption = 'Label4'
    end
    object Label5: TLabel
      Left = 320
      Top = 72
      Width = 31
      Height = 13
      Caption = 'Label5'
    end
    object LabelSum: TLabel
      Left = 384
      Top = 72
      Width = 45
      Height = 13
      Caption = 'LabelSum'
    end
    object RadioGroup1: TRadioGroup
      Left = 3
      Top = 16
      Width = 467
      Height = 41
      Caption = 'Please, Select side of block'
      Columns = 6
      ItemIndex = 0
      Items.Strings = (
        'min X (W)'
        'max X (E)'
        'min Y (S)'
        'max Y (N)'
        'min Z (B)'
        'max Z (T)')
      TabOrder = 0
      OnClick = RadioGroup1Click
    end
  end
  object GroupBox2: TGroupBox
    Left = 3
    Top = 295
    Width = 470
    Height = 82
    Caption = 
      'Input Temperature for calculating heat flux density JW,JE,JS,JN,' +
      'JB,JT.'
    TabOrder = 1
    object Label12: TLabel
      Left = 16
      Top = 24
      Width = 55
      Height = 13
      Caption = 'T min X (W)'
    end
    object Label13: TLabel
      Left = 16
      Top = 56
      Width = 55
      Height = 13
      Caption = 'T max X (E)'
    end
    object Label14: TLabel
      Left = 156
      Top = 24
      Width = 6
      Height = 13
      Caption = 'K'
    end
    object Label15: TLabel
      Left = 156
      Top = 56
      Width = 6
      Height = 13
      Caption = 'K'
    end
    object Label16: TLabel
      Left = 191
      Top = 24
      Width = 51
      Height = 13
      Caption = 'T min Y (S)'
    end
    object Label17: TLabel
      Left = 186
      Top = 56
      Width = 56
      Height = 13
      Caption = 'T max Y (N)'
    end
    object Label18: TLabel
      Left = 328
      Top = 24
      Width = 51
      Height = 13
      Caption = 'T min Z (B)'
    end
    object Label19: TLabel
      Left = 324
      Top = 51
      Width = 55
      Height = 13
      Caption = 'T max Z (T)'
    end
    object Label20: TLabel
      Left = 440
      Top = 32
      Width = 6
      Height = 13
      Caption = 'K'
    end
    object Label21: TLabel
      Left = 440
      Top = 51
      Width = 6
      Height = 13
      Caption = 'K'
    end
    object Label22: TLabel
      Left = 304
      Top = 24
      Width = 6
      Height = 13
      Caption = 'K'
    end
    object Label23: TLabel
      Left = 304
      Top = 56
      Width = 6
      Height = 13
      Caption = 'K'
    end
    object Edit7: TEdit
      Left = 77
      Top = 16
      Width = 73
      Height = 21
      TabOrder = 0
      Text = '300'
    end
    object Edit8: TEdit
      Left = 80
      Top = 48
      Width = 70
      Height = 21
      TabOrder = 1
      Text = '300'
    end
    object Edit9: TEdit
      Left = 248
      Top = 24
      Width = 50
      Height = 21
      TabOrder = 2
      Text = '300'
    end
    object Edit10: TEdit
      Left = 248
      Top = 51
      Width = 50
      Height = 21
      TabOrder = 3
      Text = '300'
    end
    object Edit11: TEdit
      Left = 385
      Top = 24
      Width = 49
      Height = 21
      TabOrder = 4
      Text = '300'
    end
    object Edit12: TEdit
      Left = 385
      Top = 51
      Width = 49
      Height = 21
      TabOrder = 5
      Text = '300'
    end
  end
  object GroupBox3: TGroupBox
    Left = 0
    Top = 383
    Width = 473
    Height = 114
    Caption = 'Result calculation density radiation heat flux'
    TabOrder = 2
    object Label24: TLabel
      Left = 16
      Top = 24
      Width = 37
      Height = 13
      Caption = 'Label24'
    end
    object Label25: TLabel
      Left = 16
      Top = 56
      Width = 37
      Height = 13
      Caption = 'Label25'
    end
    object Label26: TLabel
      Left = 189
      Top = 24
      Width = 37
      Height = 13
      Caption = 'Label26'
    end
    object Label27: TLabel
      Left = 189
      Top = 56
      Width = 37
      Height = 13
      Caption = 'Label27'
    end
    object Label28: TLabel
      Left = 341
      Top = 24
      Width = 37
      Height = 13
      Caption = 'Label28'
    end
    object Label29: TLabel
      Left = 341
      Top = 56
      Width = 37
      Height = 13
      Caption = 'Label29'
    end
    object Label30: TLabel
      Left = 136
      Top = 88
      Width = 80
      Height = 13
      Caption = 'relaxation factor'
    end
    object Button1: TButton
      Left = 362
      Top = 86
      Width = 75
      Height = 25
      Caption = 'Solve'
      TabOrder = 0
      OnClick = Button1Click
    end
    object Editrelaxfactor: TEdit
      Left = 232
      Top = 88
      Width = 121
      Height = 21
      TabOrder = 1
      Text = '0.1'
    end
    object CheckBoxinit: TCheckBox
      Left = 16
      Top = 88
      Width = 97
      Height = 17
      Caption = 'init J'
      Checked = True
      State = cbChecked
      TabOrder = 2
    end
  end
  object GroupBox4: TGroupBox
    Left = 0
    Top = 104
    Width = 473
    Height = 66
    Caption = 'Geometry dimensional Prism'
    TabOrder = 3
    object Label6: TLabel
      Left = 16
      Top = 32
      Width = 22
      Height = 13
      Caption = 'xL ='
    end
    object Label7: TLabel
      Left = 133
      Top = 32
      Width = 22
      Height = 13
      Caption = 'yL ='
    end
    object Label8: TLabel
      Left = 250
      Top = 32
      Width = 21
      Height = 13
      Caption = 'zL ='
    end
    object EditxL: TEdit
      Left = 44
      Top = 29
      Width = 66
      Height = 21
      TabOrder = 0
      Text = '0.0'
    end
    object EdityL: TEdit
      Left = 161
      Top = 29
      Width = 62
      Height = 21
      TabOrder = 1
      Text = '0.0'
    end
    object EditzL: TEdit
      Left = 277
      Top = 29
      Width = 60
      Height = 21
      TabOrder = 2
      Text = '0.0'
    end
  end
  object GroupBoxemissivity: TGroupBox
    Left = 5
    Top = 0
    Width = 467
    Height = 98
    Caption = 'emissivity individual side of block'
    TabOrder = 4
    object Label9: TLabel
      Left = 16
      Top = 24
      Width = 82
      Height = 13
      Caption = 'epsilon min X (W)'
    end
    object Label10: TLabel
      Left = 16
      Top = 56
      Width = 82
      Height = 13
      Caption = 'epsilon max X (E)'
    end
    object Label11: TLabel
      Left = 151
      Top = 24
      Width = 78
      Height = 13
      Caption = 'epsilon min Y (S)'
    end
    object Label31: TLabel
      Left = 151
      Top = 56
      Width = 83
      Height = 13
      Caption = 'epsilon max Y (N)'
    end
    object Label32: TLabel
      Left = 297
      Top = 24
      Width = 78
      Height = 13
      Caption = 'epsilon min Z (B)'
    end
    object Label33: TLabel
      Left = 293
      Top = 64
      Width = 82
      Height = 13
      Caption = 'epsilon max Z (T)'
    end
    object Edit1: TEdit
      Left = 104
      Top = 16
      Width = 41
      Height = 21
      TabOrder = 0
      Text = '0.8'
    end
    object Edit2: TEdit
      Left = 104
      Top = 58
      Width = 41
      Height = 21
      TabOrder = 1
      Text = '0.8'
    end
    object Edit3: TEdit
      Left = 235
      Top = 16
      Width = 46
      Height = 21
      TabOrder = 2
      Text = '0.8'
    end
    object Edit4: TEdit
      Left = 235
      Top = 56
      Width = 52
      Height = 21
      TabOrder = 3
      Text = '0.8'
    end
    object Edit5: TEdit
      Left = 381
      Top = 16
      Width = 67
      Height = 21
      TabOrder = 4
      Text = '0.8'
    end
    object Edit6: TEdit
      Left = 381
      Top = 56
      Width = 67
      Height = 21
      TabOrder = 5
      Text = '0.8'
    end
  end
end
