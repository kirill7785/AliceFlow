object CabinetForm: TCabinetForm
  Left = 343
  Top = 137
  Caption = 'cabinet'
  ClientHeight = 341
  ClientWidth = 613
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object PanelMain: TPanel
    Left = 8
    Top = 8
    Width = 601
    Height = 329
    Color = clMoneyGreen
    TabOrder = 0
    object GroupBoxMCS: TGroupBox
      Left = 263
      Top = 8
      Width = 329
      Height = 146
      Caption = 'Moving coordinate system'
      Color = clMoneyGreen
      ParentBackground = False
      ParentColor = False
      TabOrder = 0
      object GroupBoxangle: TGroupBox
        Left = 8
        Top = 24
        Width = 161
        Height = 117
        Caption = 'angles'
        Color = clMoneyGreen
        ParentBackground = False
        ParentColor = False
        TabOrder = 0
        object Lalfa: TLabel
          Left = 8
          Top = 20
          Width = 12
          Height = 13
          Caption = 'Alf'
        end
        object Lbeta: TLabel
          Left = 8
          Top = 52
          Width = 16
          Height = 13
          Caption = 'Bet'
        end
        object lblAlphadiap: TLabel
          Left = 96
          Top = 19
          Width = 51
          Height = 13
          Caption = 'rad [-Pi..Pi]'
        end
        object lblBetadiap: TLabel
          Left = 96
          Top = 45
          Width = 51
          Height = 13
          Caption = 'rad [-Pi..Pi]'
        end
        object LGam: TLabel
          Left = 8
          Top = 80
          Width = 22
          Height = 13
          Caption = 'Gam'
        end
        object Label1: TLabel
          Left = 96
          Top = 80
          Width = 51
          Height = 13
          Caption = 'rad [-Pi..Pi]'
        end
        object EAlf: TEdit
          Left = 41
          Top = 16
          Width = 49
          Height = 21
          TabOrder = 0
        end
        object EBet: TEdit
          Left = 41
          Top = 43
          Width = 49
          Height = 21
          TabOrder = 1
        end
        object EGam: TEdit
          Left = 40
          Top = 80
          Width = 49
          Height = 21
          TabOrder = 2
        end
      end
      object GroupBoxOrigin: TGroupBox
        Left = 184
        Top = 29
        Width = 129
        Height = 117
        Caption = 'Origin'
        Color = clMoneyGreen
        ParentBackground = False
        ParentColor = False
        TabOrder = 1
        object LabelxO: TLabel
          Left = 16
          Top = 24
          Width = 13
          Height = 13
          Caption = 'Xo'
        end
        object LabelyO: TLabel
          Left = 16
          Top = 56
          Width = 13
          Height = 13
          Caption = 'Yo'
        end
        object LabelzO: TLabel
          Left = 16
          Top = 88
          Width = 13
          Height = 13
          Caption = 'Zo'
        end
        object EditXo: TEdit
          Left = 40
          Top = 16
          Width = 81
          Height = 21
          TabOrder = 0
        end
        object EditYo: TEdit
          Left = 40
          Top = 48
          Width = 81
          Height = 21
          TabOrder = 1
        end
        object EditZo: TEdit
          Left = 43
          Top = 83
          Width = 78
          Height = 21
          TabOrder = 2
        end
      end
    end
    object GBCabinetMaterial: TGroupBox
      Left = 8
      Top = 183
      Width = 217
      Height = 129
      Caption = 'Cabinet Fluid Material'
      Color = clMoneyGreen
      ParentBackground = False
      ParentColor = False
      TabOrder = 1
      object RGCabinetMaterial: TRadioGroup
        Left = 3
        Top = 21
        Width = 121
        Height = 105
        Caption = 'Define'
        Color = clMoneyGreen
        ItemIndex = 1
        Items.Strings = (
          'Program Library'
          'User-Defined')
        ParentBackground = False
        ParentColor = False
        TabOrder = 0
        OnClick = RGCabinetMaterialClick
      end
      object BEditMaterial: TButton
        Left = 145
        Top = 21
        Width = 57
        Height = 105
        Caption = 'Edit'
        TabOrder = 1
        OnClick = BEditMaterialClick
      end
    end
    object GBOperatingTemperature: TGroupBox
      Left = 263
      Top = 263
      Width = 177
      Height = 57
      Caption = 'Operating Temperature'
      Color = clMoneyGreen
      ParentBackground = False
      ParentColor = False
      TabOrder = 2
      object LCentiGrade: TLabel
        Left = 112
        Top = 24
        Width = 14
        Height = 13
        Caption = ' '#176#1057
      end
      object EOpTemp: TEdit
        Left = 8
        Top = 24
        Width = 97
        Height = 21
        TabOrder = 0
      end
    end
    object Panel1: TPanel
      Left = 263
      Top = 160
      Width = 321
      Height = 97
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 3
      object LabelFilmCoefficient: TLabel
        Left = 116
        Top = 67
        Width = 50
        Height = 13
        Caption = 'W/(m!2*K)'
        Visible = False
      end
      object Label2: TLabel
        Left = 16
        Top = 16
        Width = 151
        Height = 13
        Caption = 'Temperature ambient conditions'
      end
      object Label3: TLabel
        Left = 16
        Top = 72
        Width = 45
        Height = 13
        Caption = 'film coeff.'
      end
      object EditFilmCoefficient: TEdit
        Left = 67
        Top = 67
        Width = 43
        Height = 21
        TabOrder = 0
        Text = '3'
        Visible = False
      end
      object ComboBoxFilmCoeff: TComboBox
        Left = 24
        Top = 40
        Width = 242
        Height = 21
        ItemIndex = 0
        TabOrder = 1
        Text = 'adiabatic wall'
        OnChange = ComboBoxFilmCoeffChange
        Items.Strings = (
          'adiabatic wall'
          'heat transfer coeff. Newton-Richman condition'
          'heat transfer coeff. Stefan-Bolcman condition'
          'heat transfer coeff. mix condition')
      end
    end
    object BCentiGrade: TButton
      Left = 471
      Top = 282
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 4
      OnClick = BCentiGradeClick
    end
    object Panel2: TPanel
      Left = 8
      Top = 8
      Width = 249
      Height = 169
      Color = clMoneyGreen
      ParentBackground = False
      TabOrder = 5
      object Label5: TLabel
        Left = 18
        Top = 16
        Width = 72
        Height = 13
        Caption = 'Geometry Type'
      end
      object ComboBoxGeometryTypeCabinet: TComboBox
        Left = 112
        Top = 8
        Width = 112
        Height = 21
        ItemIndex = 0
        TabOrder = 0
        Text = 'Prism'
        OnChange = ComboBoxGeometryTypeCabinetChange
        Items.Strings = (
          'Prism'
          'None')
      end
      object GroupBoxcabinetsize: TGroupBox
        Left = 3
        Top = 42
        Width = 238
        Height = 117
        Caption = 'cabinet size'
        Color = clMoneyGreen
        ParentBackground = False
        ParentColor = False
        TabOrder = 1
        object LxS: TLabel
          Left = 23
          Top = 27
          Width = 12
          Height = 13
          Caption = 'xS'
        end
        object LyS: TLabel
          Left = 23
          Top = 56
          Width = 12
          Height = 13
          Caption = 'yS'
        end
        object LzS: TLabel
          Left = 23
          Top = 88
          Width = 12
          Height = 13
          Caption = 'zS'
        end
        object LxE: TLabel
          Left = 127
          Top = 27
          Width = 12
          Height = 13
          Caption = 'xE'
        end
        object LyE: TLabel
          Left = 127
          Top = 56
          Width = 12
          Height = 13
          Caption = 'yE'
        end
        object LzE: TLabel
          Left = 127
          Top = 88
          Width = 12
          Height = 13
          Caption = 'zE'
        end
        object Labelunit1: TLabel
          Left = 208
          Top = 27
          Width = 16
          Height = 13
          Caption = 'mm'
        end
        object Labelunit2: TLabel
          Left = 208
          Top = 51
          Width = 16
          Height = 13
          Caption = 'mm'
        end
        object Labelunit3: TLabel
          Left = 208
          Top = 86
          Width = 16
          Height = 13
          Caption = 'mm'
        end
        object Labelunit4: TLabel
          Left = 96
          Top = 27
          Width = 16
          Height = 13
          Caption = 'mm'
        end
        object Labelunit5: TLabel
          Left = 96
          Top = 51
          Width = 16
          Height = 13
          Caption = 'mm'
        end
        object Labelunit6: TLabel
          Left = 96
          Top = 86
          Width = 16
          Height = 13
          Caption = 'mm'
        end
        object ExS: TEdit
          Left = 41
          Top = 24
          Width = 49
          Height = 21
          TabOrder = 0
        end
        object EyS: TEdit
          Left = 41
          Top = 53
          Width = 49
          Height = 21
          TabOrder = 1
        end
        object EzS: TEdit
          Left = 41
          Top = 83
          Width = 49
          Height = 21
          TabOrder = 2
        end
        object ExE: TEdit
          Left = 145
          Top = 24
          Width = 57
          Height = 21
          TabOrder = 3
        end
        object EyE: TEdit
          Left = 145
          Top = 53
          Width = 57
          Height = 21
          TabOrder = 4
        end
        object EzE: TEdit
          Left = 145
          Top = 83
          Width = 57
          Height = 21
          TabOrder = 5
        end
      end
    end
  end
end
